# A good guide: http://www.bayesian-inference.com/credible
# http://stats.stackexchange.com/questions/24681/what-is-the-decision-theoretic-justification-for-bayesian-credible-interval-proc

library(modeest)
library(coda)
s<-rnorm(10000,0,0.01)

plot_dens <- function(x, y, low=NULL, high=NULL, point=NULL, color, est_color, point_est_line_color, overdraw_interval) {
        no_estimates <- if(overdraw_interval) 1 else max(length(point), length(low))
        old_par <- par(mar=c(0.5,1,0.5,1), ylbias = 0.1)
        plot(x, y, type="l", col=rgb(0,0,0,0), yaxt = "n", xaxt = "n", ylab="", bty="n",
             ylim=c(-max(y) * 0.05 * no_estimates, max(y)), xlim = range(c(x, low, high)),xlab="")
        polygon(c(x[1] ,x, x[length(x)] ), c(0, y, 0), col="snow2", border = NA)
        #mtext("X label string", 1, line=1)
        for(i in seq_along(low)) {
                x1 <- min(which(x >= low[i]))
                x2 <- max(which(x <  high[i]))
                polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=color[i], border = NA)
                if(overdraw_interval) {
                        lines(c(low[i], high[i]), rep(-max(y) * 0.05 , 2), lwd=2*i, col = est_color[i])
                } else {
                        lines(c(low[i], high[i]), rep(-max(y) * 0.05 * i , 2), lwd=2, col = est_color[i])
                }
        }
        
        lines(x, y, type="l")
        lines(range(x), c(0,0), type="l", col="gray")
        
        for(i in seq_along(point)) {
                if(overdraw_interval) {
                        points(point[i], -max(y) * 0.05, pch=19, col = est_color[i], lwd=2 * max(1, length(low)))
                } else {
                        points(point[i], -max(y) * 0.05 * i , pch=19, col = est_color[i], lwd=2)
                        lines(c(point[i], point[i]), c(0, y[which.min(abs(x - point[i]))]), lty=2, col=point_est_line_color[i], lwd=2)
                }
        }
        par(old_par)
}


plot_and_calc_dens <- function(dist, xlim, calc_point=NULL, calc_interval=NULL, coverage = 0.95,
                               color, est_color, point_est_line_color, overdraw_interval,...) {
        calc_point <- c(calc_point)
        calc_interval <- c(calc_interval)
        n = 256
        x = seq(xlim[1], xlim[2], length.out = n)
        
        y <- get(paste0("d", dist))(x, ...)
        s <-  get(paste0("r", dist))(99999, ...)
        low <- c()
        high <- c()
        point <- c()
        coverage <- rep(coverage, length(calc_interval))
        
        for(i in seq_along(calc_interval)) {
                interval <- calc_interval[[i]](s, coverage[i])
                low[i] <- interval[1]
                high[i] <- interval[2]
        }
        for(i in seq_along(calc_point)) {
                point[i] <- calc_point[[i]](s)
        }
        
        plot_dens(x, y, low, high, point, color, est_color, point_est_line_color, overdraw_interval)
}

ddual_norm <- function(x, ...) {
        dnorm(x, -2) + dnorm(x, 2) * 0.9
}

rdual_norm <- function(n, ...) {
        c(rnorm(n / 2 * 1.1, -2), rnorm(n/2 * 0.9, 2))
}

plot_density_grid <- function(calc_point=NULL, calc_interval =NULL, coverage = 0.95, color, est_color, point_est_line_color, overdraw_interval=FALSE) {
        old_par <- par(mfrow=c(2, 3))
        plot_and_calc_dens("norm", c(-3, 3), calc_point, calc_interval, coverage = coverage, color, est_color, point_est_line_color, overdraw_interval)
        plot_and_calc_dens("exp", c(-0.1, 4), calc_point, calc_interval, rate=1, coverage = coverage, color, est_color, point_est_line_color, overdraw_interval)
        plot_and_calc_dens("beta", c(0.3, 1), calc_point, calc_interval, shape1=10, shape2=2, coverage = coverage, color, est_color, point_est_line_color, overdraw_interval)
        plot_and_calc_dens("gamma", c(0, 14), calc_point, calc_interval, shape=4, coverage = coverage, color, est_color, point_est_line_color, overdraw_interval)
        plot_and_calc_dens("lnorm", c(0, 27), calc_point, calc_interval, sd=1.2, mean=1, coverage = coverage, color, est_color, point_est_line_color, overdraw_interval)
        plot_and_calc_dens("dual_norm", c(-5, 5), calc_point, calc_interval, rate=1, coverage = coverage, color, est_color, point_est_line_color, overdraw_interval)
        par(old_par)
}

mode <- function(s) {mlv(s, method = "HSM")$M}
sd_interval <- function(s, coverage) {
        se <- sd(s) * qnorm(1 - (1 -coverage) / 2)
        c(mean(s) -  se, mean(s) + se)
}
quantile_interval <- function(s, coverage) {quantile(s, c((1 - coverage) / 2, 1 - (1 - coverage) / 2))}
hdi_interval <- function(s, coverage) {HPDinterval(mcmc(s), coverage)}

coverage <- 0.95
set.seed(123)

#plot_density_grid(mode, hdi_interval, coverage, color = "lightblue", "darkblue", "darkblue")

estimate_mode <- function(s) {
        d <- density(s)
        d$x[which.max(d$y)]
}
estimate_mode(s)

library(modeest)
mlv(s, method = "HSM")
library(coda)
HPDinterval(mcmc(s), 0.95)
plot_density_grid(median, quantile_interval, coverage, color = "lightgreen", "darkgreen", "darkgreen")

median(s)
quantile(s, c(0.025, 0.975)) # for a 95% interval.

plot_density_grid(mean, sd_interval, coverage,color =  "salmon", "darkred", "darkred")
#plot_density_grid(c(mode, median, mean), c(hdi_interval, quantile_interval, sd_interval),
#                  coverage, color = rep(rgb(0, 0, 0, 0), 3), c("darkblue", "darkgreen", "darkred"), rep(rgb(0, 0, 0, 0), 3))
plot_density_grid(median, hdi_interval, coverage, color = "orange1", "orange4", "orange4")
#plot_density_grid(c(mode), c(hdi_interval, hdi_interval), overdraw_interval = TRUE,
#                  coverage = c(0.95, 0.5), color = rep(rgb(0, 0, 0, 0), 3), c("darkred", "darkred", "darkblue"), rep(rgb(0, 0, 0, 0), 3))
limit_and_norm <- function(s, lim) {
        s <- s[ s >= lim[1] & s <= lim[2]]
        #s <- s - min(s)
        s <- s / diff(range(s))
        #s <- s + runif(1, 0, 0.9)
}

n <- 10000
d <- data.frame(dist = rep(1:6, each = n),
                s = c(
                        rnorm(n),
                        rexp(n, rate=1),
                        rbeta(n, shape1=10, shape2=2),
                        rgamma(n, shape=4),
                        rlnorm(n, sd=1.2, mean=1),
                        rdual_norm(n)
                )
)

lims <- list(c(-3, 3), c(0, 4), c(0.3, 1), c(0, 14), c(0, 27), c(-5, 5))

d <- plyr::ddply(d, "dist", function(d_sub) {
        data.frame(s = limit_and_norm(d_sub$s, lims[[ d_sub$dist[1] ]] ) )
})



qplot(x=as.factor(dist), y=s, data=d, geom="violin", fill=I("lightblue")) + theme_minimal() + theme(axis.line=element_blank(),
                                                                                                    axis.text.x=element_blank(),                                                                                            panel.grid.minor=element_blank()
                                                                                                    #plot.background=element_blank()
)
