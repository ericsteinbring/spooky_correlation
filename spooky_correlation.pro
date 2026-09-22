; load some useful functions first
forward_function quantum_correlation
forward_function bell_theorem

pro spooky_correlation

; This is a model of comprehensively tallying possible pairwise correlations of two sources such that they could be exploited to mimic a violation of a Bell test. It repeats the calculations, as a check, which reproduce the closed-form versions intended for the original journal paper. Various flavours of those calculations can be tried, to verify they are correct.

; settings
M = 0.87 ; 1. ; 0.99 ; 0.89 ; 0.88300 ; 0.87 ; 0.866 ; sqrt(3./4.) ; 0.86 ; 0.84 ; 0.75
;R = 0.66 ; 1. ; 0.68 ; 0.67 ; 0.667 ; 0.66667 ; 2./3. ; 0.661 ; 0.66 ; 0.5 ; 0.25
R_discriminant = 0.25 ; 1. ; 0.5 ; 0.25 ; 0.
fold_over = 'yes' ; 'no' ; 'yes' ; whether to "fold" the observed distribution, looking for absolute differences, as was reported in the original paper
sampling = 2e6 ; 1e6 ; fine-scaling of the sampling of the residual distribution
resolution = 1000 ; floor(sqrt(2e6)) ; 1000 ; 500 ; 100 ; resolution of calculations
trials = 50000 ; 30000 ; 20000 ; 10000 ; 1000
reduced = 0.2 ; 0.01 ; 0.1 ; 0.2 ; 0.5 ; the fraction of detections
check_original = 'no' ; 'no' ; 'yes' ; just show the original equation, and stop
just_distribution = 'no' ; 'no' ; 'yes' ; only show the residual distribution, and stop

; display setup
loadct, 0, /silent ; greyscale, for models
;loadct, 39, /silent ; colour, for data
!p.background=255
device, set_font='helvetica bold', /tt_font
!p.font=1
!except=0
charsize=3.5 ; 3. ; 2.5 ; 2.

; set up a display
window, 0, xsize = 1000, ysize = 800, title="Bell's theorem versus expected quantum-mechanical correlation"
wset, 0
!p.multi = [0, 1, 1]

; plotting up the Bell-theorem and expected quantum-mechanical correlation functions
theta = 180.*findgen(resolution)/float(resolution) ; degrees
p = double(bell_theorem(theta))
q = double(quantum_correlation(theta))
plot, theta, p, xstyle=1, xrange=[-5., 185.], xtitle='Discriminator angle (deg)', ystyle=1, yrange=[-1., 2.], ytitle='Correlation of measurements', charsize=charsize, color=0
oplot, theta, q, thick=2, color=0
oplot, theta, p, color=0
; limits
oplot, [-5., 185.], [0., 0.], linestyle=1, color=0
;oplot, [-5., 185.], [0.5, 0.5], linestyle=1, color=0
oplot, [-5., 185.], [mean(p), mean(p)], linestyle=2, color=0
oplot, [-5., 185.], [median(p), median(p)], linestyle=3, color=0
; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[-5., 185.], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[-1., 2.], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[-1., 2.], ytickformat='(A1)'

; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_bell_versus_quantum_correlation.jpg', screen_output, quality=100, true=1

; report the maximal differences
D = sqrt(2.)*2. ; Tirelson bound
print, 'Maximal difference (Tirelson bound): ' , D
D = max(abs(p))+max(abs(q))
print, 'Maximal difference (both p and q max): ', D
D = max(abs(p - q))
print, 'Maximal difference (max of abs(p-q)): ', D
D = max(sqrt(p^2. - q^2.))
print, 'Maximal difference (geometric max.): ', D

; starting over once again from scratch, and possibly rewritten as converted to using an input of M only
phi = 180.*findgen(resolution)/float(resolution) ; degrees
p = double(bell_theorem(phi))
q = double(quantum_correlation(phi))
;R = 2./3.*(1./float(resolution))*total(M*sqrt(1. + (4./3.)*p^2. - (4./3.)*q^2.))
;R = 2./3.*(1./float(resolution))*total(M*sqrt(1. + (1./M^2.)*p^2. - (1./M^2.)*q^2.)) ; converted
;R = 2.*(1./!pi)*(1./float(resolution))*total(M*sqrt(1. + (4./3.)*p^2. - (4./3.)*q^2.)) ; changed
;R = 2.*(1./!pi)*(1./float(resolution))*total(M*sqrt(1. + (1./M^2.)*p^2. - (1./M^2.)*q^2.)) ; corrected
R = 2.*(1./!pi)*(1./float(resolution))*total(M*sqrt(1. + (4./3.)*p^2. - (4./3.)*q^2.)) ; changed back to original
print, 'M = ', M
print, 'R = ', R
;C = (3./4.)*(1./float(resolution))*total(sqrt(1. + p^2. - q^2.)) ; special case of D where R=0
;C = (3./4.)*(1./float(resolution))*total(sqrt(1. + (4./3.)*p^2. - (4./3.)*q^2.)) ; special case
C = (3./4.)^2.*(1./float(resolution))*total(sqrt(1. + (4./3.)*p^2. - (4./3.)*q^2.)) ; change back to original
;C = M^4.*(1./float(resolution))*total(sqrt(1. + p^2. - q^2.)) ; converted, not universal
;C = M^4.*(1./float(resolution))*total(sqrt(1. + (1./M^2.)*p^2. - (1./M^2.)*q^2.)) ; corrected
;M_0 = C - sqrt(3./4.)
M_0 = C - M ; converted
print, 'M, R, C, M_0;', M, R, C, M_0
;D = sqrt(M^2. + p^2. - q^2.)
D = (3./2.)*(M*sqrt(1. + (4./3.)*p^2. - (4./3.)*q^2.) - R) ; changed back to original
;D = M + (p - q)
;D = abs(M + p - q)
;D = M + abs(p - q)
D_R = (3./2.)*(D - R)
;D_R = !pi/2.*(D - R) ; changed
D_0 = 0.5*(abs(D_R) - R)
D_p = 0.5*(D_R - R)
;D_0_p = abs(median(D_0))/4.
;D_l = 3.*(abs(D_R) - (4./3.)*2.*R)
;D_l = 3.*(abs(D_R) - (1./M^2.)*2.*R) ; converted
D_l = 3.*(abs(D) - (4./3.)*2.*R) ; changed back to original
D_u = abs(D_R) + R
D_t = (D_l + D_u)/2.
D_p = 0.5*(D_R - R)
if fold_over eq 'yes' then D_t(where(abs(D_R) lt 2.*R)) = D_0(where(abs(D_R) lt 2.*R))
P = 0.5 + 0.5*D_t

; zeropoint
;D_0_p = 2.*M_0
D_0_p = abs(median(D_0))/4.
print, 'D_0_p: ', D_0_p

;; hold onto these
;P_hold = P
;P_0_hold = P_0
;P_0_p_hold = P_0_p

; calculating the probabilities directly
;P = 0.5 + 0.5*D
;P = 0.5 + 0.5*D_0
P = 0.5 + 0.5*D_0 - D_0_p
;if fold_over eq 'yes' then P(where(abs(D_R) ge 2.*R)) = 0.5 + 2.*D_0(where(abs(D_R) ge 2.*R)) - (3./4.)*R
if fold_over eq 'yes' then P(where(abs(D_R) ge 2.*R)) = 0.5 + 2.*D_0(where(abs(D_R) ge 2.*R)) - (M^2.)*R ; converted
P_0 = 0.5 + 0.5*D_0
P_0_p = 0.5 + 0.5*D_p

; inflections
phi_0_i = min(phi(where(P ge R)))
phi_0_o = max(phi(where(P ge R)))
print, 'Inflections (degrees): ', phi_0_i, phi_0_o, phi_0_o-phi_0_i
D_phi_0 = (7./4.)*R
phi_0_i = min(phi(where(P_0_p le min(P))))
phi_0_o = max(phi(where(P_0_p le min(P))))
P_phi_0 = 0.5 - (D_phi_0 - 1.)
print, 'D_phi_0, P_phi_0: ', D_phi_0, P_phi_0 
print, 'Cusps (degrees): ', phi_0_i, phi_0_o, phi_0_o - phi_0_i

; crossing
phi_0 = max(phi(where(P ge median(P))))
print, 'Crossing median (degrees): ', phi_0
phi_0 = max(phi(where(P ge mean(P))))
print, 'Crossing mean (degrees): ', phi_0

; hold onto these
P_hold = P
P_0_hold = P_0
P_0_p_hold = P_0_p

; set up a display
window, 1, xsize = 1000, ysize = 800, title='Distribution of residuals'
wset, 1
!p.multi = [0, 1, 1]

; generate a randomized distribution of residuals
;R_trial = 0.333*randomn(seed, sampling)+0.5 ; normal distibution of trials
;R_trial = 0.333*randomn(seed, sampling)+0.67 ; normal distibution of trials
;R_trial = 0.5*randomn(seed, sampling)+0.5 ; normal distibution of trials
;R_trial = 0.5*randomu(seed, sampling)+0.5 ; uniform distribution of trials
R_trial = R*randomn(seed, sampling) + 0.5 ; normal distibution of trials
;R_trial(where(R_trial gt 1.)) = !values.f_nan
;R_trial(where(R_trial ge 0.5)) = !values.f_nan
;R_trial(where(R_trial lt 0.)) = !values.f_nan
;R_trial = R_trial(sort(R_trial))
;binsize = (max(R_trial)-min(R_trial))/float(sampling)
;binsize = (max(R_trial)-min(R_trial))/(float(sampling)/float(trials))
binsize = (max(R_trial)-min(R_trial))/float(resolution)
distribution_R = histogram(R_trial, min=min(R_trial), max=max(R_trial), binsize=binsize, /nan)
distribution_R = smooth(distribution_R, 100, /edge_wrap)
distribution_R = distribution_R/float(max(distribution_R))
;distribution_R = distribution_R/float(total(distribution_R))
;distribution_R = 0.5*distribution_R/float(max(distribution_R))
bins_R = binsize*findgen(n_elements(distribution_R)) + min(R_trial)
distribution_R = (0.25/distribution_R(min(where(bins_R ge 0.))))*distribution_R ; setting it to have 25% probabability at zero
distribution_R = (R_discriminant/distribution_R(min(where(bins_R ge 0.))))*distribution_R ; setting it to have probabability set to the discrimant value at zero
print, 'Min, max, mean, median of R: ', min(R_trial), max(R_trial), mean(R_trial), median(R_trial), R_trial(min(where(distribution_R eq max(distribution_R))))

; plot this distribution up
plot, bins_R, distribution_R, xstyle=1, xrange=[min(R_trial), max(R_trial)], xtitle='Possible residual', ystyle=1, yrange=[0., 1.1*max(distribution_R)], ytitle='Relative probability', charsize=charsize, color=0
; shading
for i = 0, n_elements(bins_R) - 1 do begin
 oplot, [bins_R(i), bins_R(i)], [0., distribution_R(i)], thick=2, color=200
 if bins_R(i) gt 0. and bins_R(i) lt 1. then begin 
  oplot, [bins_R(i), bins_R(i)], [0., distribution_R(i)], thick=2, color=150
 endif
endfor
print, 'Relative probability at 0: ', distribution_R(min(where(bins_R ge 0.)))
print, 'Relative probability at 1: ', distribution_R(min(where(bins_R ge 1.)))
; limits
;oplot, [0.5, 0.5], [0., 1.1*max(distribution_R)], thick=5, color=100
oplot, [0.5, 0.5], [0., max(distribution_R)], thick=5, color=100
oplot, [mean(R_trial), mean(R_trial)], [0., 1.1*max(distribution_R)], linestyle=2, color=0
oplot, [median(R_trial), median(R_trial)], [0., 1.1*max(distribution_R)], linestyle=3, color=0
oplot, [bins_R(min(where(distribution_R eq max(distribution_R)))), bins_R(min(where(distribution_R eq max(distribution_R))))], [0., 1.1*max(distribution_R)], linestyle=1, color=0
; replot
oplot, bins_R, distribution_R, color=0
;; some more limits
;oplot, [0., 0.], [0., 1.1*max(distribution_R)], color=0
;oplot, [1., 1.], [0., 1.1*max(distribution_R)], color=0
oplot, [min(R_trial), max(R_trial)], [R_discriminant, R_discriminant], linestyle=1, color=0

; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[min(R_trial), max(R_trial)], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[min(R_trial), max(R_trial)], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[0., 1.1*max(distribution_R)], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[0., 1.1*max(distribution_R)], ytickformat='(A1)'

; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_residual_distribution.jpg', screen_output, quality=100, true=1
if just_distribution eq 'yes' then stop

; set up the display
window, 2, xsize = 1000, ysize = 800, title='Probability calculations'
wset, 2
!p.multi = [0, 1, 1]

; plot this
plot, phi, 0.*P, xstyle=1, xrange=[-5., 185.], xtitle='Separation angle (deg)', ystyle=1, yrange=[0., 1.05], ytitle='Probability of correlated flux', charsize=charsize, color=0
oplot, phi, P, color=0
;oplot, phi, P_0, color=0
;oplot, phi, P_0_p, color=0
;oplot, phi, p, color=0
;oplot, phi, q, color=0
;oplot, phi, P, thick=3, color=100
if check_original eq 'yes' then stop

; and perform a simulation of a randomized distribution of residuals 
P_trial = fltarr(resolution)
R_trial = congrid(R_trial, resolution) ; gridding the sampling of the distribution down to the resolution of the calculations
tally = 0
outsides = 0
detections = 0
for i = 0, resolution - 1 do begin
 p = double(bell_theorem(phi))
 q = double(quantum_correlation(phi))
; D = abs(M + p - q)
 D = (3./2.)*(sqrt(M^2. + p^2. - q^2.) - R_trial(i))
; D = (!pi/2.)*(sqrt(M^2. + p^2. - q^2.) - R_trial(i)) ; changed
 D_0 = 0.5*(abs(D) - mean(abs(D)))
; P = 0.5 + 0.5*D_0
 P = 0.5 + 0.5*D_0 - D_0_p
; if fold_over eq 'yes' then P(where(abs(D) ge 2.*mean(abs(D)))) = 0.5 + 2.*D_0(where(abs(D) ge 2.*mean(abs(D)))) - mean(abs(D))
; if fold_over eq 'yes' then P(where(abs(D) ge 2.*mean(abs(D)))) = 0.5 + 2.*D_0(where(abs(D) ge 2.*mean(abs(D)))) - (3./4.)*mean(abs(D))
; if fold_over eq 'yes' then P(where(abs(D) ge 2.*mean(abs(D)))) = 0.5 + 2.*D_0(where(abs(D) ge 2.*mean(abs(D)))) - (!pi/4.)*mean(abs(D)) ; changed
 if fold_over eq 'yes' then P(where(abs(D) ge 2.*mean(abs(D)))) = 0.5 + 2.*D_0(where(abs(D) ge 2.*mean(abs(D)))) - (M^2.)*mean(abs(D)) ; converted
 fraction = randomu(seed, resolution) ; a random distribution from zero to unity
 reduced_fraction = fraction
 fraction(where(fraction lt 1.-float(trials)/float(resolution)^2.)) = 0. ; throw out some
 fraction(where(fraction ge 1.-float(trials)/float(resolution)^2.)) = 1. ; keep some
 reduced_fraction(where(fraction lt 1.-float(reduced*trials)/float(resolution)^2.)) = 0. ; throw out some
 reduced_fraction(where(fraction ge 1.-float(reduced*trials)/float(resolution)^2.)) = 1. ; keep some
 tally = tally + total(fraction)
 outside = total(fraction(where(fraction*P(where(phi gt 90.)) gt 0.75)))
 outsides = outsides + outside
 detected = total(fraction(where(fraction*P gt 0.75)))
 detections = detections + detected
 plotsym, 0, 0.5, /fill, color=200
 oplot, phi, fraction*P, psym=8, color=0
 if (R_trial(i) gt 0. and R_trial(i) lt 1.) then begin
  plotsym, 0, 0.5, /fill, color=150
  oplot, phi, fraction*P, psym=8, color=0
;  if R_trial(i) eq 0.5 then begin
  if (R_trial(i) gt 0.45 and R_trial(i) lt 0.55) then begin
   plotsym, 0, 0.5, /fill, color=100
   oplot, phi, fraction*P, psym=8, color=0
  endif
 endif
; plotsym, 0, 0.5, /fill, color=0
; oplot, phi, reduced_fraction*P, psym=8, color=0
; P_trial = P_trial + fraction*P
 P_trial = P_trial + P
endfor
P_trial = P_trial/float(resolution) ;+ D_0_p ; /float(trials)
;P_trial(where(P_trial lt 0.)) = !values.f_nan
P_trial(where(phi eq max(phi))) = P_trial(where(phi eq max(phi))-1)
oplot, phi, P_trial, thick=3, color=100
print, 'Trials, total detections, those outside horizon, and their fractions:'
print, tally, detections, outsides, float(detections)/float(tally), float(outsides)/float(tally)
print, 'Giving cases inside/outside of:', tally*float(detections)/float(tally), tally*float(outsides)/float(tally)
print, 'And those with 10 samples each: ', tally*float(detections)/float(tally)/20., tally*float(outsides)/float(tally)/20.
print, 'And those with 20 samples each of: ', tally*float(detections)/float(tally)/40., tally*float(outsides)/float(tally)/40.

; limits
;oplot, phi, replicate(M, resolution), color=0
;oplot, phi, replicate(R, resolution), color=0
;oplot, phi, replicate(0.75, resolution), linestyle=3, color=0
oplot, phi, replicate(0.75, resolution), linestyle=1, color=0
;oplot, phi, replicate(0.75, resolution), color=0
oplot, phi, replicate(0.5, resolution), linestyle=3, color=0
oplot, phi, replicate(0.25, resolution), linestyle=1, color=0
;oplot, phi, replicate(0.25, resolution), color=0
oplot, phi, replicate(mean(P), resolution), linestyle=2, color=0

; replot
;oplot, phi, P_0_p_hold, color=0
;oplot, phi, P_0_hold, color=0
oplot, phi, P_hold, color=0

; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[-5., 185.], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[0., 1.05], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[0., 1.05], ytickformat='(A1)'

; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_probability_calculation.jpg', screen_output, quality=100, true=1

; report
print, 'Max P, P_0, P_trial:', max(P_hold), max(P_0_hold), max(P_trial, /nan)
print, 'Mean P, P_0, P_trial:', mean(P_hold), mean(P_0_hold), mean(P_trial, /nan)
print, 'Median P, P_0, P_trial:', median(P_hold), median(P_0_hold), median(P_trial)
;print, 'Antipode P, P_0, P_trial:', P_hold(where(phi eq max(phi))), P_0_hold(where(phi eq max(phi))), P_trial(where(phi eq max(phi, /nan)))

; proportion
proportion = (1.-max(P_hold))/mean(P_hold)
print, 'Proportion: ', proportion

; signal-to-noise ratio
alpha_sky = 0.250
beta_sky  = 0.400
gamma_sky = 0.375
epsilon = 0.250
zeta = 0.02 ; 0.02 ; 0.015 ; 0.0125 ; 0.01
omega = 1.000
print, 'Samples: ', 1
; ideal case
snr_limit_space = omega/zeta ; proportion*omega/zeta
snr_long = proportion*omega/((1./sqrt(float(1.)))*(alpha_sky + beta_sky*gamma_sky^2.+epsilon) + zeta)/snr_limit_space
print, 'Idealized S/N (space, long): ', snr_limit_space, snr_long
; normalized to space
snr_limit_space = proportion*omega/zeta
snr_long = omega/((1./sqrt(float(1.)))*(alpha_sky + beta_sky*gamma_sky^2.+epsilon) + zeta)/snr_limit_space
print, 'Normalized S/N (space, long): ', snr_limit_space, snr_long
samples = trials
print, 'Samples: ', samples
; ideal case
snr_limit_space = omega/zeta ; proportion*omega/zeta
snr_long = proportion*omega/((1./sqrt(float(samples)))*(alpha_sky + beta_sky*gamma_sky^2.+epsilon) + zeta)/snr_limit_space
print, 'Idealized S/N (space, long): ', snr_limit_space, snr_long
; normalized to space
snr_limit_space = proportion*omega/zeta
snr_long = omega/((1./sqrt(float(samples)))*(alpha_sky + beta_sky*gamma_sky^2.+epsilon) + zeta)/snr_limit_space
print, 'Normalized S/N (space, long): ', snr_limit_space, snr_long
samples = 2.e6 ; sampling
print, 'Samples: ', samples
; ideal case
snr_limit_space = omega/zeta ; proportion*omega/zeta
snr_long = proportion*omega/((1./sqrt(float(samples)))*(alpha_sky + beta_sky*gamma_sky^2.+epsilon) + zeta)/snr_limit_space
print, 'Idealized S/N (space, long): ', snr_limit_space, snr_long
; normalized to space
snr_limit_space = proportion*omega/zeta
snr_long = omega/((1./sqrt(float(samples)))*(alpha_sky + beta_sky*gamma_sky^2.+epsilon) + zeta)/snr_limit_space
print, 'Normalized S/N (space, long): ', snr_limit_space, snr_long

end

; quantum correlation
function quantum_correlation, theta
; Returns the "classical" two-point quantum correlation of point sources. Theta is in degrees.
;theta = theta-180. ; test of mirroring
;return, -1.*(-1.*cos(!pi*theta/180.)) ; test of mirroring
return, -1.*cos(!pi*theta/180.)
end

; Bell theorem
function bell_theorem, theta
; Returns the two-point quantum correlation of point sources, as in Bell's Theorem. Theta is in degrees.
;theta = theta-180. ; test of mirroring
;return, -1.*abs(cos(2.*!pi*theta/180.) - cos(!pi*theta/180.)) + cos(!pi*theta/180.) ; test of mirroring
return, abs(cos(2.*!pi*theta/180.) - cos(!pi*theta/180.)) + cos(!pi*theta/180.)
end
