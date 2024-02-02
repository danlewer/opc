# power estimate for overdose prevention centres in Sandwell and Cardiff

# ----------------------
# sources of assumptions
# ----------------------

# hospital admissions
# https://digital.nhs.uk/data-and-information/publications/statistical/statistics-on-public-health/2023/data-tables
# 22.7/100,000 admissions with primary diagnosis of poisoning by drug misuse
# -> 56m / 100,000 * 22.7 = 12,712 admissions per year
# -> 34 per day

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013357/
# 4% of calls to ambulance call centre in Wales 2007/8 coded as overdose 

# https://aace.org.uk/uk-ambulance-service/national-ambulance-data/
# approx. 30,000 calls to ambulance control centres per day (England I think)
# -> 1,200 overdose-related calls per day

# drug-related deaths
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/drugmisusedeathsbylocalauthority
# 4907 deaths registered in England & Wales 2022
# 16 in Sandwell
# 22 in Cardiff

# estimated ambulance call outs per LA - scaled by drug-related deaths
# sandwell: 16 / 4907 * 56/(56+3) * 1200
# -> 4/day or 113/month
# Cardiff: 22 / 4907 * 56/(56+3) * 1200
# -> 5/day or 155/month
# 35000 per month in England
# 36500 per month in England & Wales

# -----------------------------
# function to do one simulation
# -----------------------------

# before - number of periods before the intervention
# after - number of periods after the intervention
# effect - rate ratio
# basei - baseline count of events per period before the intervention in the intervention area
# basec - count of events per period in the control area
# assume no secular trend

sim <- function(before = 36, after = 12, effect = 0.9, basei = 113, basec = 35000) {
  cts <- rpois(before + after, basec)
  its <- rpois(before + after, rep(c(basei, basei * effect), times = c(before, after)))
  outcome <- c(cts, its)
  time <- rep(1:(before + after), times = 2)
  intervention <- rep(0:1, each = before + after)
  before_after <- rep(rep(c('before', 'after'), times = c(before, after)), 2)
  m <- glm(outcome ~ time + intervention*before_after, family = 'poisson')
  summary(m)$coef[5,4]
}

# -----------------------------------------------------
# function to do multiple simulations and esimate power
# -----------------------------------------------------

power <- function(B = 5000, ...) {
  p <- sapply(seq_len(B), function (x) {
    if (x %% 100 == 0) print (x)
    sim (...)
  })
  mean(p < 0.05)
}

# --------------
# estimate power
# --------------

effects <- seq(0.75, 1, 0.01)
power_estimates <- sapply(effects, function (x) power(effect = x))

# find 80% power

log_effects <- log(effects)
m <- glm(power_estimates ~ log_effects, family = 'binomial')
d <- data.frame(log_effects = seq(-0.3, 0, 0.001))
d$power <- predict(m, newdata = d, type = 'response')
log_effect0.8 <- d[which.min(abs(d$power - 0.8)),'log_effects']
exp(log_effect0.8)

# plot

xs <- c(0.75, 0.8, 0.85, 0.9, 0.95, 1)

png('opc_power.png', height = 6, width = 6, units = 'in', res = 300)

plot(1, type = 'n', xlim = c(-0.25, 0), ylim = c(0, 1), axes = F, xlab = 'Risk ratio', ylab = 'Power')
points(log_effects, power_estimates)
with(d[d$log_effects >= -0.25,], lines(log_effects, power))
segments(-0.25, 0.8, x1 = log_effect0.8, lty = 2, col = 'red')
segments(log_effect0.8, 0, y1 = 0.8, lty = 2, col = 'red')
rect(-0.25, 0, 0, 1)
axis(2, 0:5/5, paste0(0:5 * 20, '%'), pos = -0.25, las = 2)
axis(1, log(xs), xs, pos = 0)
axis(1, log_effect0.8, labels = round(exp(log_effect0.8), 2), tick = F, pos = 0, col.axis = 'red')

dev.off()
