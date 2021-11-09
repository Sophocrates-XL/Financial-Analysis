library(stable);

# This namespace includes some wrapper functions that either makes the use
# of functions in the stable package more convenient, or performs batch
# analysis for model selection.
Stable = list(
    
    # Calculates the loglikelihood for a stablereg model.
    Loglik = function(model) {
        y = model$y;
        output = 0;
        for (i in 1:length(y)) {
            output = output + log(dstable(y[i], tail = model$tail[i],
                skew = model$skew[i], disp = model$disp[i], loc = model$loc[i]));
        }
        return(output);
    },
    
    # Calculates the AIC for a stablereg model.
    AIC = function(model) {
        model.loglik = Stable$Loglik(model);
        n.est = length(model$coefficients);
        return(2 * n.est - 2 * model.loglik);
    },
    
    # Calculates the BIC for a stablereg model.
    BIC = function(model) {
        y = model$y;
        model.loglik = Stable$Loglik(model);
        n.sample = length(y);
        n.est = length(model$coefficients);
        return(n.est * log(n.sample) - 2 * model.loglik);
    },
    
    # Calculates the quantiles for a ONE-PERIOD stablereg model
    # over time duration t.
    # The default is to calculate the deciles.
    Quantiles = function(model, t, probabilities = seq(from = 0.9, to = 0.1, by = -0.1)) {
        tail = model$tail[1];
        skew = model$skew[1];
        disp = model$disp[1] * sqrt(t);
        loc = model$loc[1] * t;
        log.quantiles = qstable(probabilities, tail = tail, skew = skew, disp = disp, loc = loc);
        quantiles = exp(log.quantiles) - 1;
        return(data.frame(
            probability = probabilities,
            return = quantiles
        ));
    },
    
    # Calculates the VaR at p for a ONE-PERIOD stablereg model.
    VaR = function(model, t, p) {
        tail = model$tail[1];
        skew = model$skew[1];
        disp = model$disp[1] * sqrt(t);
        loc = model$loc[1] * t;
        return(exp(qstable(p, tail = tail, skew = skew, disp = disp, loc = loc)) - 1);
    },
    
    # Provides a list of all possible ONE-PERIOD stablereg models with tail, skew, and disp
    # constrained or unconstrained.
    # Initial values for optimization of tail, skew and disp are set according to
    # the setup of the original publisher of the stable package, because of
    # computational reasons.
    # Returns a list that includes a list of all models and a data frame of model statistics.
    Stablereg.List = function(sample) {
        models = list(
            center.normal = stablereg(y = sample, tail = tail_g(1.9999999), otail = F,
                skew = 0, oskew = F, disp = ~ 1, idisp = -3, loc = 0, oloc = F),
            sym.center.stable = stablereg(y = sample, tail = ~ 1, itail = 1.9999,
                skew = 0, oskew = F, disp = ~ 1, idisp = -3, loc = 0, oloc = F),
            center.stable = stablereg(y = sample, tail = ~ 1, itail = 1.9999, skew = ~ 1,
                iskew = 0, disp = ~ 1, idisp = -3, loc = 0, oloc = F),
            normal = stablereg(y = sample, tail = tail_g(1.9999999), otail = F, skew = 0,
                iskew = 0, disp = ~ 1, idisp = -3, loc = ~ 1, iloc = 0),
            sym.stable = stablereg(y = sample, tail = ~ 1, itail = 1.9999, skew = 0,
                oskew = F, disp = ~ 1, idisp = -3, loc = ~ 1, iloc = 0),
            stable = stablereg(y = sample, tail = ~ 1, itail = 1.9999, skew = ~ 1,
                iskew = 0, disp = ~ 1, idisp = -3, loc = ~ 1, iloc = 0)
        );
        tails = rep(0, length(models));
        skews = rep(0, length(models));
        disps = rep(0, length(models));
        locs = rep(0, length(models));
        logliks = rep(0, length(models));
        AICs = rep(0, length(models));
        BICs = rep(0, length(models));
        for (i in 1:length(models)) {
            tails[i] = models[[i]]$tail[1];
            skews[i] = models[[i]]$skew[1];
            disps[i] = models[[i]]$disp[1];
            locs[i] = models[[i]]$loc[1];
            logliks[i] = Stable$Loglik(models[[i]]);
            AICs[i] = Stable$AIC(models[[i]]);
            BICs[i] = Stable$BIC(models[[i]]);
        }
        model.names = names(models);
        print(sprintf("Model with minimum AIC: %s.",
            model.names[which(AICs == min(AICs))]));
        print(sprintf("Model with minimum BIC: %s.",
            model.names[which(BICs == min(BICs))]));
        return(list(
            models = models,
            result.df = data.frame(
                model = model.names,
                tail = tails,
                skew = skews,
                disp = disps,
                loc = locs,
                loglik = logliks,
                AIC = AICs,
                BIC = BICs
            )
        ));
    }
            
);



# Stable model for coal ETF.
coal = read.csv("Coal ETF.csv", header = T);
coal.close = coal$C;
plot(coal.close);
ln.coal.close = log(coal.close);
ln.coal.diff = ln.coal.close[2:length(ln.coal.close)] - ln.coal.close[1:(length(ln.coal.close) - 1)];
plot(ln.coal.diff, cex = 0.5, col = "blue");
lines(rep(0, length(ln.coal.diff)), lty = "dashed");

coal.model.set = Stable$Stablereg.List(ln.coal.diff);
coal.model.set$result.df;
coal.preferred.model = coal.model.set$models$normal;

ln.coal.diff.1 = ln.coal.diff[1:240];
ln.coal.diff.2 = ln.coal.diff[241:length(ln.coal.diff)];
plot(ln.coal.diff, cex = 0.5, col = c(rep("blue", 240), rep("red", length(ln.coal.diff) - 240)),
     xlab = "t", ylab = "ln(X[t + 1]/X[t])",
     main = "Coal ETF: Historical performance");
lines(rep(0, length(ln.coal.diff)), lty = "dashed");
coal.1.model.set = Stable$Stablereg.List(ln.coal.diff.1);
coal.2.model.set = Stable$Stablereg.List(ln.coal.diff.2);
coal.1.model.set$result.df;
coal.2.model.set$result.df;
coal.1.preferred.model = coal.1.model.set$models$normal;
coal.2.preferred.model = coal.2.model.set$models$normal;
coal.combined.loglik = Stable$Loglik(coal.1.preferred.model) +
    Stable$Loglik(coal.2.preferred.model);
coal.combined.AIC = 2 * 4 - 2 * coal.combined.loglik;
coal.combined.BIC = 4 * log(length(ln.coal.diff)) - 2 * coal.combined.loglik;
coal.combined.loglik;
coal.initial.loglik = Stable$Loglik(coal.preferred.model);
coal.initial.loglik;
2 * (coal.combined.loglik - coal.initial.loglik);
pchisq(2 * (coal.combined.loglik - coal.initial.loglik), df = 2, lower.tail = F);
coal.combined.AIC;
Stable$AIC(coal.preferred.model);
coal.combined.BIC;
Stable$BIC(coal.preferred.model);



# Stable model for S&P 500 ETF.
sp500 = read.csv("S&P ETF.csv", header = T);
sp.close = sp500$C;
plot(sp.close);
ln.sp.close = log(sp.close);
ln.sp.diff = ln.sp.close[2:length(ln.sp.close)] - ln.sp.close[1:(length(ln.sp.close) - 1)];
plot(ln.sp.diff, xlab = "t", ylab = "ln(X[t + 1]/X[t])",
     main = "S&P500 ETF: Historical performance", cex = 0.5, col = "blue");
lines(rep(0, length(ln.sp.diff)), lty = "dashed");
sp.model.set = Stable$Stablereg.List(ln.sp.diff);
sp.model.set$result.df;
Stable$Quantiles(sp.model.set$models$stable, t = 250);
Stable$Quantiles(sp.model.set$models$normal, t = 250);
VaR.t = c(30, 60, 120, 240);
stable.VaRs.5 = rep(0, length(VaR.t));
normal.VaRs.5 = rep(0, length(VaR.t));
stable.VaRs.1 = rep(0, length(VaR.t));
normal.VaRs.1 = rep(0, length(VaR.t));
for (i in 1:length(VaR.t)) {
    stable.VaRs.5[i] = Stable$VaR(sp.model.set$models$stable, t = VaR.t[i], p = 0.05);
    normal.VaRs.5[i] = Stable$VaR(sp.model.set$models$normal, t = VaR.t[i], p = 0.05);
    stable.VaRs.1[i] = Stable$VaR(sp.model.set$models$stable, t = VaR.t[i], p = 0.01);
    normal.VaRs.1[i] = Stable$VaR(sp.model.set$models$normal, t = VaR.t[i], p = 0.01);
}
stable.VaRs.5;
normal.VaRs.5;
stable.VaRs.1;
normal.VaRs.1;
