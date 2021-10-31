BSM = list(
    
    C = function(S, K, t, r, s) {
        if (K == 0) {
            return(S);
        } else {
            d1 = (log(S / K) + (r + s^2 / 2) * t) / (s * sqrt(t));
            d2 = d1 - s * sqrt(t);
            return(pnorm(d1) * S - pnorm(d2) * K * exp(-r * t));
        }
    }

);

Valuation = list(
    
    # Construct a template of the OPM valuation model which stores the company's
    # capital structure to facilitate reusability over different parameter values.
    Construct = function(thresholds, allocation.matrix, FDS) {
        
        return(list(
            
            PPS = function(equity.value, time.left, risk.free.rate, volatility) {
                n.tranche = length(thresholds);
                n.share.type = nrow(allocation.matrix);
                call.values = rep(0, n.tranche);
                for (i in 1:length(call.values)) {
                    call.values[i] = BSM$C(S = equity.value, K = thresholds[i],
                        t = time.left, r = risk.free.rate, s = volatility);
                };
                V = diag(call.values);
                V.star = diag(c(call.values[2:n.tranche], 0));
                A = allocation.matrix;
                if (sum(is.na(FDS)) == 0) {
                    FDS.matrix = diag(FDS);
                    return(solve(FDS.matrix) %*% A %*% (V - V.star) %*% rep(1, n.tranche));
                } else {
                    aggregate.values = A %*% (V - V.star) %*% rep(1, n.tranche);
                    return(ifelse(is.na(FDS), NA, aggregate.values / FDS));
                }
            },
            
            PPS.Delta = function(equity.value, time.left, risk.free.rate, volatility) {
                n.tranche = length(thresholds);
                n.share.type = nrow(allocation.matrix);
                S = equity.value;
                t = time.left;
                r = risk.free.rate;
                s = volatility;
                d1.vec = log(S / thresholds) / (s * sqrt(t)) + r * sqrt(t) / s +
                    s * sqrt(t) / 2;
                DELTA = diag(pnorm(d1.vec));
                DELTA.star = diag(c(pnorm(d1.vec[2:n.tranche]), 0));
                A = allocation.matrix;
                if (sum(is.na(FDS)) == 0) {
                    FDS.matrix = diag(FDS);
                    return(solve(FDS.matrix) %*% A %*% (DELTA - DELTA.star) %*% rep(1, n.tranche));
                } else {
                    deriv.aggregate.values = A %*% (DELTA - DELTA.star) %*% rep(1, n.tranche);
                    return(ifelse(is.na(FDS), NA, deriv.aggregate.values / FDS));
                }
            },
            
            PPS.Vega = function(equity.value, time.left, risk.free.rate, volatility) {
                n.tranche = length(thresholds);
                n.share.type = nrow(allocation.matrix);
                S = equity.value;
                t = time.left;
                r = risk.free.rate;
                s = volatility;
                d1.vec = log(S / thresholds) / (s * sqrt(t)) + r * sqrt(t) / s +
                    s * sqrt(t) / 2;
                NU = diag(S * sqrt(t) * dnorm(d1.vec));
                NU.star = diag(S * sqrt(t) * c(dnorm(d1.vec[2:n.tranche]), 0));
                A = allocation.matrix;
                if (sum(is.na(FDS)) == 0) {
                    FDS.matrix = diag(FDS);
                    return(solve(FDS.matrix) %*% A %*% (NU - NU.star) %*% rep(1, n.tranche));
                } else {
                    deriv.aggregate.values = A %*% (NU - NU.star) %*% rep(1, n.tranche);
                    return(ifelse(is.na(FDS), NA, deriv.aggregate.values / FDS));
                }
            }
            
        ));

    }

);



# Case 1 concerns the hypothetical case of Cotopaxi Tech.
FDS = c(5000000, NA, 15000000, 2000000, 1000000, 10000000, 3000000);
thresholds = c(0, 7500000, 13500000, 33000000, 38250000, 42650000, 58750000, 91750000);
allocation.matrix = rbind(
    c(1, 0, 0, 0.714, 0.625, 0.217, 0.152, 0.139),
    c(0, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0.652, 0.455, 0.417),
    c(0, 0, 0, 0.286, 0.25, 0.087, 0.061, 0.056),
    c(0, 0, 0, 0, 0.125, 0.043, 0.03, 0.028),
    c(0, 0, 0, 0, 0, 0, 0.303, 0.278),
    c(0, 0, 0, 0, 0, 0, 0, 0.083)
);
n.tranche = ncol(allocation.matrix);
n.share.type = nrow(allocation.matrix);
case.1.model = Valuation$Construct(thresholds = thresholds,
    allocation.matrix = allocation.matrix, FDS = FDS);
case.1.model$PPS(equity.value = 40000000, time.left = 3, risk.free.rate = 0.02, volatility = 0.8);



# Sensitivity analysis for valuation.
colors = c("red", "pink", "green", "blue", "purple", "black");
legends = c("Series B preferred shares", "Series A preferred shares", "Common shares",
    "Options", "Warrants I", "Warrants II");

valuations = seq(from = 20000000, to = 60000000, by = 100000);
n.valuation.trial = length(valuations);

PPS.matrix = matrix(data = NA, nrow = n.valuation.trial, ncol = n.share.type);
for (i in 1:n.valuation.trial) {
    PPS.matrix[i,] = case.1.model$PPS(equity.value = valuations[i], time.left = 3,
                                      risk.free.rate = 0.02, volatility = 0.8);
}
PPS.plotable.matrix = PPS.matrix[,c(1, 3:n.share.type)];
colors = c("red", "pink", "green", "blue", "purple", "black");
plot(PPS.plotable.matrix[,1] ~ valuations, ylim = c(0, max(PPS.plotable.matrix)),
     type = "l", ylab = "PPS", xlab = "equity valuation",
     main = "Cotopaxi Tech: PPS over equity valuation");
for (i in 1:ncol(PPS.plotable.matrix)) {
    lines(PPS.plotable.matrix[,i] ~ valuations, col = colors[i]);
}
lines(rep(0, n.valuation.trial) ~ valuations, lty = "dashed");
legend(x = "topleft", legend = legends, y.intersp = 0.25, col = colors,
       lty = "solid", bty = "n");



# Sensitivity analysis for volatility.
volatilities = seq(from = 0.4, to = 1.2, by = 0.001);
n.volatility.trial = length(volatilities);

PPS.matrix = matrix(data = NA, nrow = n.volatility.trial, ncol = n.share.type);
for (i in 1:n.volatility.trial) {
    PPS.matrix[i,] = case.1.model$PPS(equity.value = 40000000, time.left = 3,
        risk.free.rate = 0.02, volatility = volatilities[i]);
}
PPS.plotable.matrix = PPS.matrix[,c(1, 3:n.share.type)];
plot(PPS.plotable.matrix[,1] ~ volatilities, ylim = c(0, max(PPS.plotable.matrix)),
     type = "l", ylab = "PPS", xlab = "volatility",
     main = "Cotopaxi Tech: PPS over volatility");
for (i in 1:ncol(PPS.plotable.matrix)) {
    lines(PPS.plotable.matrix[,i] ~ volatilities, col = colors[i]);
}
lines(rep(0, n.volatility.trial) ~ volatilities, lty = "dashed");
legend(x = "topleft", legend = legends, y.intersp = 0.25, col = colors,
    lty = "solid", bty = "n");

volatilities = seq(from = 0.1, to = 2, by = 0.001);
n.volatility.trial = length(volatilities);
PPS.vega.matrix = matrix(data = NA, nrow = n.volatility.trial, ncol = n.share.type);
for (i in 1:n.volatility.trial) {
    PPS.vega.matrix[i,] = case.1.model$PPS.Vega(equity.value = 40000000, time.left = 3,
        risk.free.rate = 0.02, volatility = volatilities[i]);
}
PPS.vega.plotable.matrix = PPS.vega.matrix[,c(1, 3:n.share.type)];
plot(PPS.vega.plotable.matrix[,1] ~ volatilities,
     ylim = c(min(PPS.vega.plotable.matrix), max(PPS.vega.plotable.matrix)),
     type = "l", ylab = "PPS vega", xlab = "volatility",
     main = "Cotopaxi Tech: PPS vega over volatility");
for (i in 1:ncol(PPS.vega.plotable.matrix)) {
    lines(PPS.vega.plotable.matrix[,i] ~ volatilities, col = colors[i]);
}
lines(rep(0, n.volatility.trial) ~ volatilities, lty = "dashed");
legend(x = "topleft", legend = legends, y.intersp = 0.25, col = colors,
       lty = "solid", bty = "n");



# Case 2 concerns the hypothetical case of ComplexCo.
thresholds = c(0, 2500, 5500, 11500, 23000);
allocation.matrix = rbind(
    c(0.4, 0, 0.25, 0.217, 0.196),
    c(0.6, 0, 0, 0.13, 0.118),
    c(0, 1, 0.75, 0.652, 0.588),
    c(0, 0, 0, 0, 0.098)
);
FDS = c(500, 300, 1500, 250);
n.tranche = length(thresholds);
n.share.type = nrow(allocation.matrix);
equity.value = 17500;
time.left = 4;
volatility = 0.35;
risk.free.rate = 0.015;
case.2.model = Valuation$Construct(thresholds = thresholds,
    allocation.matrix = allocation.matrix, FDS = FDS);
case.2.model$PPS(equity.value = equity.value, time.left = time.left,
    volatility = volatility, risk.free.rate = risk.free.rate);



# Sensitivity analysis for valuation.
colors = c("red", "green", "blue", "black");
legends = c("Class A preferred shares", "Class B preferred shares", "Common shares",
    "Warrants");

valuations = seq(from = 10000, to = 25000, by = 50);
n.valuation.trial = length(valuations);

PPS.matrix = matrix(data = NA, nrow = n.valuation.trial, ncol = n.share.type);
for (i in 1:n.valuation.trial) {
    PPS.matrix[i,] = case.2.model$PPS(equity.value = valuations[i], time.left = time.left,
        risk.free.rate = risk.free.rate, volatility = volatility);
}
PPS.plotable.matrix = PPS.matrix;
plot(PPS.plotable.matrix[,1] ~ valuations, ylim = c(0, max(PPS.plotable.matrix)),
     type = "l", ylab = "PPS", xlab = "equity valuation",
     main = "ComplexCo: PPS over equity valuation");
for (i in 1:ncol(PPS.plotable.matrix)) {
    lines(PPS.plotable.matrix[,i] ~ valuations, col = colors[i]);
}
lines(rep(0, n.valuation.trial) ~ valuations, lty = "dashed");
legend(x = "topleft", legend = legends, y.intersp = 0.25, col = colors,
       lty = "solid", bty = "n");


# Sensitivity analysis for volatilities.
colors = c("red", "green", "blue", "black");
legends = c("Class A preferred shares", "Class B preferred shares", "Common shares",
            "Warrants");

volatilities = seq(from = 0.1, to = 0.6, by = 0.001);
n.volatility.trial = length(volatilities);

PPS.matrix = matrix(data = NA, nrow = n.volatility.trial, ncol = n.share.type);
for (i in 1:n.volatility.trial) {
    PPS.matrix[i,] = case.2.model$PPS(equity.value = equity.value, time.left = time.left,
                                      risk.free.rate = risk.free.rate, volatility = volatilities[i]);
}
PPS.plotable.matrix = PPS.matrix;
plot(PPS.plotable.matrix[,1] ~ volatilities, ylim = c(0, max(PPS.plotable.matrix)),
     type = "l", ylab = "PPS", xlab = "volatility",
     main = "ComplexCo: PPS over volatility");
for (i in 1:ncol(PPS.plotable.matrix)) {
    lines(PPS.plotable.matrix[,i] ~ volatilities, col = colors[i]);
}
lines(rep(0, n.volatility.trial) ~ volatilities, lty = "dashed");
legend(x = "topleft", legend = legends, y.intersp = 0.25, col = colors,
       lty = "solid", bty = "n");

volatilities = seq(from = 0.1, to = 2, by = 0.001);
n.volatility.trial = length(volatilities);

PPS.vega.matrix = matrix(data = NA, nrow = n.volatility.trial, ncol = n.share.type);
for (i in 1:n.volatility.trial) {
    PPS.vega.matrix[i,] = case.2.model$PPS.Vega(equity.value = equity.value, time.left = time.left,
        risk.free.rate = risk.free.rate, volatility = volatilities[i]);
}
PPS.vega.plotable.matrix = PPS.vega.matrix;
plot(PPS.vega.plotable.matrix[,1] ~ volatilities,
     ylim = c(min(PPS.vega.plotable.matrix), max(PPS.vega.plotable.matrix)),
     type = "l", ylab = "PPS vega", xlab = "volatility",
     main = "ComplexCo: PPS vega over volatility");
for (i in 1:ncol(PPS.vega.plotable.matrix)) {
    lines(PPS.vega.plotable.matrix[,i] ~ volatilities, col = colors[i]);
}
lines(rep(0, n.volatility.trial) ~ volatilities, lty = "dashed");
legend(x = "topleft", legend = legends, y.intersp = 0.25, col = colors,
       lty = "solid", bty = "n");

# Sensitivity analysis for a specific position.
w = c(0.5, 0.25, 0.25, 0);
volatilities = seq(from = 0.1, to = 0.6, by = 0.001);
n.volatility.trial = length(volatilities);
mean.share.prices = rep(0, n.volatility.trial);
for (i in 1:n.volatility.trial) {
    p = case.2.model$PPS(equity.value = equity.value, volatility = volatilities[i],
        risk.free.rate = risk.free.rate, time.left = time.left);
    mean.share.prices[i] = as.numeric(t(w) %*% p);
}
plot(PPS.plotable.matrix[,1] ~ volatilities, ylim = c(0, max(PPS.plotable.matrix)),
     type = "l", ylab = "PPS", xlab = "volatility",
     main = "ComplexCo: Mean price of share combination over individual share types");
for (i in 1:ncol(PPS.plotable.matrix)) {
    lines(PPS.plotable.matrix[,i] ~ volatilities, col = colors[i]);
}
lines(mean.share.prices ~ volatilities, type = "l", col = "purple", lwd = 3);
lines(rep(0, n.volatility.trial) ~ volatilities, lty = "dashed");
legend(x = "topleft", legend = c(legends, "50:25:25:0"), y.intersp = 0.25,
    col = c(colors, "purple"), lty = "solid", bty = "n", lwd = c(rep(1, length(legends)), 3));


w = c(0.5, 0.25, 0.25, 0);
volatilities = seq(from = 0.1, to = 2, by = 0.001);
n.volatility.trial = length(volatilities);
mean.share.price.vegas = rep(0, n.volatility.trial);
for (i in 1:n.volatility.trial) {
    p = case.2.model$PPS.Vega(equity.value = equity.value, volatility = volatilities[i],
        risk.free.rate = risk.free.rate, time.left = time.left);
    mean.share.price.vegas[i] = as.numeric(t(w) %*% p);
}
plot(PPS.vega.plotable.matrix[,1] ~ volatilities,
     ylim = c(min(PPS.vega.plotable.matrix), max(PPS.vega.plotable.matrix)),
     type = "l", ylab = "PPS vega", xlab = "volatility",
     main = "ComplexCo: Mean price vega of share combination over individual share types");
for (i in 1:ncol(PPS.vega.plotable.matrix)) {
    lines(PPS.vega.plotable.matrix[,i] ~ volatilities, col = colors[i]);
}
lines(mean.share.price.vegas ~ volatilities, type = "l", col = "purple", lwd = 3);
lines(rep(0, n.volatility.trial) ~ volatilities, lty = "dashed");
legend(x = "topleft", legend = c(legends, "50:25:25:0"), y.intersp = 0.25,
       col = c(colors, "purple"), lty = "solid", bty = "n", lwd = c(rep(1, length(legends)), 3));
