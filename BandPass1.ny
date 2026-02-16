;nyquist plug-in
;version 4
;type process
;name "Bandpass Filter (IIR/FFT, Phase, Gain)"
;action "Bandpass filtering..."

;control fc        "Center Frequency (Hz)"    real   "" 1000 0.001 22000
;control bw        "Bandwidth (Hz)"           real   "" 100  0.001 22000
;control qmode     "Q Mode"                   choice "Proportionate-Q (BW Hz),Constant-Q (BW as Q)" 0
;control design    "Filter Design"            choice "Cascaded IIR (bandpass2),True Spectral (FFT)" 0
;control phasemode "Phase Mode"               choice "Natural,Zero-Phase (FFT),Linear-Phase (FFT)" 0
;control stages    "Steepness (Stages)"       int    "" 4    1 16
;control mix       "Wet Mix (%)"              int    "" 100  0 100
;control gainmode  "Output Gain Mode"         choice "Match Envelope (RMS+Peak),Auto Comp (sqrt(Q)),User Gain (dB),Peak Normalize" 0
;control gaindb    "User Gain (dB)"           real   "" 0   -24 24

;; -------- input selection (plugin + prompt) --------
(setq sig
  (cond
    ((and (boundp 's) (soundp s)) s)
    ((and (boundp '*track*) (soundp *track*)) *track*)
    (t (error "No input sound found (expected s or *track*)."))))

;; -------- parameter safety clamps --------
(setq sr (snd-srate sig))
(setq nyq (/ sr 2.0))

(setq fc (max 0.001 (min fc (- nyq 0.001))))

(setq bw (max 0.001 bw))
(setq maxbw (* 2.0 (min fc (- nyq fc))))
(setq bw (min bw (max 0.001 maxbw)))

;; Q interpretation:
;; - Proportionate-Q: bw is bandwidth in Hz (Q varies with fc)
;; - Constant-Q: bw parameter is treated as Q directly
(setq q
  (if (= qmode 0)
      (max 0.001 (/ fc bw))
      (max 0.001 bw)))

;; compute effective BW (Hz) from Q for later (FFT edges + display/debug)
(setq effbw (max 0.001 (/ fc q)))

(setq low (max 0.001 (- fc (/ effbw 2.0))))
(setq high (min (- nyq 0.001) (+ fc (/ effbw 2.0))))

;; Optional debug print (shows in Nyquist debug output):
;; (print (list "Effective Q" q "Effective BW(Hz)" effbw))

;; -------- IIR design (cascaded bandpass2) --------
(defun do-iir (x)
  (let ((y x))
    (dotimes (k stages)
      (setq y (bandpass2 y fc q)))
    y))

;; -------- FFT design (true spectral bandpass) --------
;;
;; Uses sa-init (documented) for streaming frames, wraps it with an iterator
;; that applies a hard bin mask, and reconstructs using snd-ifft.
;;
(defun do-fft (x)
  (let* (;; choose bin width based on effective BW for "surgical" control:
         ;; finer bins when BW is narrow; clamp to keep it reasonable
         (binw (max 0.5 (min 200.0 (/ effbw 12.0))))
         (fftdur (/ 1.0 binw))
         ;; clamp window duration so it doesn't explode on ultra-narrow BW
         (fftdur (max 0.01 (min 0.50 fftdur)))
         (skipper (/ fftdur 2.0))
         (sa (sa-init :fft-dur fftdur :skip-period skipper :window :hann :input (snd-copy x)))
         (binhz (sa-get-bin-width sa))
         (fftsize (sa-get-fft-size sa))
         (halfbins (+ (/ fftsize 2) 1))
         ;; build a class that responds to :next with masked frames
         (bpclass (send class :new '(sa binhz low high phasemode fftsize halfbins))))
    ;; initializer
    (send bpclass :answer :isnew '(sa0 binhz0 low0 high0 ph0 fft0 hb0)
      '((setq sa sa0)
        (setq binhz binhz0)
        (setq low low0)
        (setq high high0)
        (setq phasemode ph0)
        (setq fftsize fft0)
        (setq halfbins hb0)
        self))
    ;; :next method
    (send bpclass :answer :next '()
      '((let ((fr (sa-next sa)))
          (if (null fr)
              nil
              (progn
                ;; zero bins outside [low, high]
                (let ((i 0))
                  (while (< i halfbins)
                    (let ((freq (* i binhz)))
                      (when (or (< freq low) (> freq high))
                        (cond
                          ((= i 0)
                           (setf (aref fr 0) 0.0))
                          ((= i (/ fftsize 2))
                           (setf (aref fr (- fftsize 1)) 0.0))
                          (t
                           (setf (aref fr (- (* 2 i) 1)) 0.0)
                           (setf (aref fr (* 2 i)) 0.0)))))
                    (setq i (+ i 1))))
                ;; phase modes (FFT only):
                ;; 0 Natural: keep original phase
                ;; 1 Zero-Phase: set imag=0, real=magnitude
                ;; 2 Linear-Phase: same as zero-phase, but add N/2 delay via (-1)^i
                (cond
                  ((= phasemode 1)
                   (let ((i 1))
                     ;; bins 1..N/2-1 have real+imag
                     (while (< i (/ fftsize 2))
                       (let* ((r (aref fr (- (* 2 i) 1)))
                              (im (aref fr (* 2 i)))
                              (mag (sqrt (+ (* r r) (* im im)))))
                         (setf (aref fr (- (* 2 i) 1)) mag)
                         (setf (aref fr (* 2 i)) 0.0))
                       (setq i (+ i 1)))))
                  ((= phasemode 2)
                   (let ((i 1))
                     (while (< i (/ fftsize 2))
                       (let* ((r (aref fr (- (* 2 i) 1)))
                              (im (aref fr (* 2 i)))
                              (mag (sqrt (+ (* r r) (* im im))))
                              ;; delay N/2 samples = multiply bin i by exp(-j*pi*i) = (-1)^i
                              (sgn (if (= (rem i 2) 0) 1.0 -1.0)))
                         (setf (aref fr (- (* 2 i) 1)) (* sgn mag))
                         (setf (aref fr (* 2 i)) 0.0))
                       (setq i (+ i 1)))
                     ;; Nyquist bin (real only) also gets (-1)^(N/2)
                     (let ((sgn (if (= (rem (/ fftsize 2) 2) 0) 1.0 -1.0)))
                       (setf (aref fr (- fftsize 1)) (* sgn (aref fr (- fftsize 1))))))))
                fr)))))
    ;; build iterator instance, reconstruct
    (let* ((it (send bpclass :new sa binhz low high phasemode fftsize halfbins))
           (skip-samps (sa-get-fft-skip-size sa))
           (y (snd-ifft (local-to-global 0) sr it skip-samps nil)))
      y)))

;; If user asks for FFT-only phase modes but design is IIR, force FFT internally
(setq use-fft (or (= design 1) (and (= design 0) (/= phasemode 0))))

;; -------- apply chosen design --------
(setq wet (if use-fft (do-fft sig) (do-iir sig)))

;; -------- output gain modes --------
(cond
  ;; 0) Match envelope (RMS follower) + peak correction (accuracy improvement)
  ((= gainmode 0)
    (let* ((rate 200.0)
           (env-in  (rms sig rate))
           (env-wet (rms wet rate))
           (scale1  (mult env-in (recip (sum env-wet 1.0e-6)))))
      (setq wet (mult wet scale1))
      ;; peak correction pass (matches overall peak, helps “accuracy”)
      (let* ((pin (peak sig ny:all))
             (pwet (peak wet ny:all))
             (s2 (/ pin (+ pwet 1.0e-9))))
        (setq wet (mult s2 wet)))))

  ;; 1) Auto compensation (audible on white noise): sqrt(Q), capped
  ((= gainmode 1)
    (setq wet (mult (min 4.0 (sqrt q)) wet)))

  ;; 2) User gain (0 dB is effectively “none”)
  ((= gainmode 2)
    (setq wet (mult (db-to-linear gaindb) wet)))

  ;; 3) Peak normalize (to ~0.99 to reduce clipping risk)
  (t
    (let* ((p (peak wet ny:all))
           (g (/ 0.99 (+ p 1.0e-9))))
      (setq wet (mult g wet)))))

;; -------- wet/dry mix (sound arithmetic via mult/sum) --------
(setq m (/ mix 100.0))
(setq out (sum (mult m wet)
               (mult (- 1.0 m) sig)))

;; -------- preview-safe: force start time = 0 --------
(cue out)
