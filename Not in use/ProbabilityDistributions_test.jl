using Statistics
using Distributions
using StatsPlots
plotly()

a = [1.0 50.0; 4.0 40; 5 20; 2 50]
b = log.(a)

# Univariat distributions
norm_I = Normal(mean(a[:, 1]), std(a[:, 1]))
norm_P = Normal(mean(a[:, 2]), std(a[:, 2]))
log_I = LogNormal(mean(b[:, 1]), std(b[:, 1]))
log_P = LogNormal(mean(b[:, 2]), std(b[:, 2]))

# Plot distributions
StatsPlots.plot(norm_I, label = "I Norm", fill = (0, 0.5, :blue))
StatsPlots.plot!(log_I, label = "I Log", fill = (0, 0.5, :orange))
StatsPlots.plot!(norm_P, label = "P Norm", fill = (0, 0.5, :blue))
StatsPlots.plot!(log_P, label = "P Log", fill = (0, 0.5, :orange))

# Sample and plot inflow samples
normI_s = rand(norm_I, 10000)
logI_s = rand(log_I, 10000)
StatsPlots.histogram(normI_s, normalize = :pdf, label = "norm", color = :blue, alpha = 0.3)
StatsPlots.histogram!(logI_s, normalize = :pdf, label = "log", color = :purple, alpha = 0.3)

# Sample and plot price samples
normP_s = rand(norm_P, 10000)
logP_s = rand(log_P, 10000)
StatsPlots.histogram(normP_s, normalize = :pdf, label = "norm", color = :blue, alpha = 0.3)
StatsPlots.histogram!(logP_s, normalize = :pdf, label = "log", color = :purple, alpha = 0.3)

# Multivariat distributions
mvnorm1 = fit(MvNormal, a')
mvnorm2 = MvNormal(vec(mean(a, dims = 1)), cov(a, dims = 1))
mvlog2 = MvLogNormal(vec(mean(b, dims = 1)), cov(b, dims = 1))

# Sample from distributions
mvNormF_s = rand(mvnorm1, 10000)
mvNorm_s = rand(mvnorm2, 10000)
mvLog_s = rand(mvlog2, 10000)

# Plot inflow samples
StatsPlots.histogram(
    mvNormF_s[1, :],
    normalize = :pdf,
    bins = 50,
    label = "samples,fitted normal",
    color = :blue,
    alpha = 0.3,
)
StatsPlots.histogram(
    mvNorm_s[1, :],
    normalize = :pdf,
    bins = 50,
    label = "samples, normal dist.",
    color = :green,
    alpha = 0.3,
)
StatsPlots.histogram!(
    mvLog_s[1, :],
    normalize = :pdf,
    bins = 50,
    label = "samples, LogNormal",
    color = :purple,
    alpha = 0.3,
)
StatsPlots.histogram!(
    a[:, 1],
    normalize = :pdf,
    label = "input data",
    color = :yellow,
    alpha = 0.3,
)

# Plot price samples
StatsPlots.histogram(
    mvNormF_s[2, :],
    normalize = :pdf,
    label = "samples,fitted normal",
    color = :blue,
    alpha = 0.3,
)
StatsPlots.histogram(
    mvNorm_s[2, :],
    normalize = :pdf,
    bins = 50,
    label = "samples, normal dist.",
    color = :green,
    alpha = 0.3,
)
StatsPlots.histogram!(
    mvLog_s[2, :],
    bins = 50,
    normalize = :pdf,
    label = "samples, LogNormal",
    color = :purple,
    alpha = 0.3,
)
StatsPlots.histogram!(
    a[:, 2],
    normalize = :pdf,
    label = "input data",
    color = :yellow,
    alpha = 0.3,
)
