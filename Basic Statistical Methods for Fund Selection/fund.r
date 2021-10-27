Plot.Manager = function(manager) {
    plot(manager$fund ~ manager$average);
};



# Order = 0: Null model.
# Order = 1: Linear model.
# Order = 2: Quadratic model with centralized quadratic term.
WLS.Manager = function(manager, order = 2) {
    fund = manager$fund;
    average = manager$average;
    duration = manager$end - manager$start;
    y.bar = log(fund / 100 + 1) / duration;
    x.bar = log(average / 100 + 1) / duration;
    if (order == 0) {
        model = lm(y.bar ~ 1, weights = duration);
    } else if (order == 1) {
        model = lm(y.bar ~ x.bar, weights = duration);
    } else if (order == 2) {
        model = lm(y.bar ~ x.bar +
            I((x.bar - mean(x.bar))^2), weights = duration);
    }
    return(model);
};



WLS.Select = function(manager) {
    duration = manager$end - manager$start;
    fund = manager$fund;
    average = manager$average;
    y.bar = log(fund / 100 + 1) / duration;
    x.bar = log(average / 100 + 1) / duration;
    null.model = WLS.Manager(manager, order = 0);
    linear.model = WLS.Manager(manager, order = 1);
    quadratic.model = WLS.Manager(manager, order = 2);
    x.display = seq(from = min(x.bar), to = max(x.bar),
        by = (max(x.bar) - min(x.bar)) / 1000);
    y.display.null = rep(coef(null.model)[1], length(x.display));
    y.display.linear = predict(linear.model, newdata = list(x.bar = x.display),
        weights = rep(1, length(x.display)));
    y.display.quadratic = predict(quadratic.model, newdata = list(x.bar = x.display),
        weights = rep(1, length(x.display)));
    y.display.max = max(y.bar, y.display.linear, y.display.quadratic);
    y.display.min = min(y.bar, y.display.linear, y.display.quadratic);
    plot(y.bar ~ x.bar, ylab = "ln(fund)", xlab = "ln(market)",
        ylim = c(y.display.min, y.display.max),
        main = sprintf("%s: Model selection with WLS estimators", manager$name));
    lines(y.display.null ~ x.display, col = "black", lwd = 2);
    lines(y.display.linear ~ x.display, col = "blue", lwd = 2);
    lines(y.display.quadratic ~ x.display, col = "purple", lwd = 2);
    legend(x = "topleft", legend = c("Null", "Linear", "Quadratic"),
        col = c("black", "blue", "purple"), lty = rep("solid", 3), lwd = 2,
        x.intersp = 0.25, y.intersp = 0.25, bty = "n");
    print(sprintf("Null model: BIC = %f.", BIC(null.model)));
    print(sprintf("Linear model: BIC = %f.", BIC(linear.model)));
    print(sprintf("Quadratic model: BIC = %f.", BIC(quadratic.model)));
    if (BIC(null.model) <= min(BIC(linear.model), BIC(null.model))) {
        print("Preferred model: null.");
    } else if (BIC(linear.model) <= BIC(quadratic.model)) {
        print("Preferred model: linear.");
    } else {
        print("Preferred model: quadratic.");
    }
    return(list(
        null.model = null.model,
        linear.model = linear.model,
        quadratic.model = quadratic.model
    ));
};



# Does NOT cover the null model, since the null model's predictions are
# independent of the covariates.
WLS.Predict = function(manager, order = 2, level = 0.8) {
    duration = manager$end - manager$start;
    fund = manager$fund;
    average = manager$average;
    y.bar = log(fund / 100 + 1) / duration;
    x.bar = log(average / 100 + 1) / duration;
    model = WLS.Manager(manager, order = order);
    x.display = seq(from = min(x.bar), to = max(x.bar),
        by = (max(x.bar) - min(x.bar)) / 1000);
    y.display.CI = predict(model, newdata = list(x.bar = x.display),
        weights = rep(1, length(x.display)), level = level, interval = "confidence");
    y.display = y.display.CI[,1];
    y.display.CI.upper = y.display.CI[,3];
    y.display.CI.lower = y.display.CI[,2];
    y.display.PI = predict(model, newdata = list(x.bar = x.display),
        weights = rep(1, length(x.display)), level = level,interval = "prediction");
    y.display.PI.upper = y.display.PI[,3];
    y.display.PI.lower = y.display.PI[,2];
    y.display.max = max(y.bar, y.display.PI.upper);
    y.display.min = min(y.bar, y.display.PI.lower);
    plot(y.bar ~ x.bar, ylab = "ln(fund)", xlab = "ln(market)",
         ylim = c(y.display.min, y.display.max),
         main = sprintf("%s: Annual %d%% CI and PI", manager$name, as.integer(level * 100)));
    lines(rep(0, length(x.display)) ~ x.display, col = "black", lty = "dashed");
    lines(y.display ~ x.display, col = "black", lwd = 2);
    lines(y.display.CI.lower ~ x.display, col = "blue", lty = "dashed");
    lines(y.display.CI.upper ~ x.display, col = "blue", lty = "dashed");
    lines(y.display.PI.lower ~ x.display, col = "purple", lty = "dashed");
    lines(y.display.PI.upper ~ x.display, col = "purple", lty = "dashed");
    legend(x = "topleft", legend = c("Predicted", "CI", "PI"), col = c("black", "blue", "purple"),
        lty = c("solid", "dashed", "dashed"), lwd = c(2, 1, 1), x.intersp = 0.25, y.intersp = 0.25,
        bty = "n");
};



GLS.Manager = function(manager, order = 2, cov.factor = 1, plot.fit = T) {
    fund = manager$fund;
    average = manager$average;
    start = manager$start;
    end = manager$end;
    duration = end - start;
    n.data = length(fund);
    y.bar = log(fund / 100 + 1) / duration;
    x.bar = log(average / 100 + 1) / duration;
    V = matrix(data = NA, nrow = n.data, ncol = n.data);
    for (i in 1:nrow(V)) {
        for (j in 1:ncol(V)) {
            if (i == j) {
                V[i, j] = 1 / duration[i];
            } else {
                V[i, j] = max(0, min(end[i], end[j]) - max(start[i], start[j])) /
                    (duration[i] * duration[j]) * cov.factor;
            }
        }
    };
    if (order == 0) {
        X = rep(1, n.data);
    } else if (order == 1) {
        X = cbind(rep(1, n.data), x.bar);
    } else if (order == 2) {
        X = cbind(rep(1, n.data), x.bar, (x.bar - mean(x.bar))^2);
    }
    beta.GLS = solve(t(X) %*% solve(V) %*% X) %*% t(X) %*% solve(V) %*% y.bar;
    if (plot.fit) {
        plot(y.bar ~ x.bar);
        x.display = seq(from = min(x.bar), to = max(x.bar),
            by = (max(x.bar) - min(x.bar)) / 1000);
        if (order == 0) {
            y.display = beta.GLS[1] + 0 * x.display;
        } else if (order == 1) {
            y.display = beta.GLS[1] + beta.GLS[2] * x.display;
        } else if (order == 2) {
            y.display = beta.GLS[1] + beta.GLS[2] * x.display + beta.GLS[3] *
                (x.display - mean(x.bar))^2;
        }
        lines(y.display ~ x.display, col = "red");
    }
    return(list(
        type = ifelse(order == 0, "null", ifelse(order == 1, "linear", "quadratic")),
        x.bar = x.bar,
        y.bar = y.bar,
        beta.GLS = beta.GLS,
        y.bar.fitted = X %*% beta.GLS
    ));
};



GLS.Plot.Family = function(manager, order = 2, moderated.only = F) {
    cov.factors = seq(from = 0, to = 1, by = 0.2);
    colors = c("black", "purple", "pink", "blue", "green", "red");
    duration = manager$end - manager$start;
    y.bar = log(manager$fund / 100 + 1) / duration;
    x.bar = log(manager$average / 100 + 1) / duration;
    plot(y.bar ~ x.bar);
    x.display = seq(from = min(x.bar), to = max(x.bar),
        by = (max(x.bar) - min(x.bar)) / 1000);
    y.display.matrix = matrix(data = NA, nrow = length(cov.factors), ncol = length(x.display));
    if (moderated.only) {
        cov.factor.i = c(2:(length(cov.factors) - 1), 1);
    } else {
        cov.factor.i = c(2:(length(cov.factors) - 1), 1, length(cov.factors));
    }
    for (i in cov.factor.i) {
        gls.model = GLS.Manager(manager, order = order,
            cov.factor = cov.factors[i], plot.fit = F);
        beta.GLS = gls.model$beta.GLS;
        if (order == 0) {
            y.display.matrix[i,] = beta.GLS[1] + 0 * x.display;
        } else if (order == 1) {
            y.display.matrix[i,] = beta.GLS[1] + beta.GLS[2] * x.display;
        } else if (order == 2) {
            y.display.matrix[i,] = beta.GLS[1] + beta.GLS[2] * x.display + beta.GLS[3] *
                (x.display - mean(x.bar))^2;
        }
    }
    plot(y.bar ~ x.bar, ylab = "ln(fund)", xlab = "ln(market)",
         ylim = c(min(y.display.matrix[cov.factor.i,], y.bar),
            max(y.display.matrix[cov.factor.i,], y.bar)),
         main = sprintf("%s: Sensitivity for GLS estimators", manager$name));
    for (i in cov.factor.i) {
        lines(y.display.matrix[i,] ~ x.display, col = colors[i],
            lty = ifelse(i == 1 | i == length(cov.factors), "solid", "dashed"),
            lwd = ifelse(i == 1 | i == length(cov.factors), 2, 1));
    };
    legend.texts = rep("", length(cov.factor.i));
    legend.colors = rep("", length(cov.factor.i));
    for (i in sort(cov.factor.i)) {
        legend.texts[i] = paste0("z = ", cov.factors[i]);
        legend.colors[i] = colors[i];
    }
    legend(x = "topleft", legend = legend.texts, col = legend.colors,
        x.intersp = 0.25, y.intersp = 0.25, lwd = c(2, 1, 1, 1, 1, 2)[sort(cov.factor.i)],
        lty = c("solid", rep("dashed", 4), "solid")[sort(cov.factor.i)], bty = "n");
};
    


now = 2021.81;
liwei = list(
    name = "Li Wei",
    average = c(-1.47, 0.62, 18.96, 94.47, 18.07, 75.81, 128.42, 55.09, 60.87, 37.42, 249.41),
    fund = c(-0.46, -0.3, 21.31, 153.67, -0.5, 172.6, 232.9, 46.1, 23.47, 3.09, 559.4),
    start = c(2021.56, 2021.06, 2020.69, 2019.41, 2017.61, 2016.14, 2016, 2014.49, 2014.21,
        2013.69, 2011.63),
    end = c(now, now, now, now, 2020.11, 2021.28, now, 2019.17, 2019.39, 2015.13, now)
);
length(liwei$fund);
Plot.Manager(liwei);
liwei.selection = WLS.Select(liwei);
summary(liwei.selection$linear.model);
GLS.Plot.Family(liwei, order = 1);
WLS.Predict(liwei, order = 1);

# This section validates if the lm results are consistent with the formulations
# in the write-up.
n = length(liwei$fund);
duration = liwei$end - liwei$start;
X = cbind(rep(1, n), log(liwei$average / 100 + 1) / duration);
y = log(liwei$fund / 100 + 1) / duration;
W = diag(duration);
beta.WLS = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y;
# Consistent with parameter estimates in lm.
beta.WLS;
SSE = as.numeric(t(y) %*% W %*% y - t(y) %*% W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y);
s2 = SSE / (n - 2);
s = sqrt(s2);
# Consistent with residual standard error in lm.
s;
# Consistent with vcov in lm.
vcov.beta = s2 * solve(t(X) %*% W %*% X);
vcov.beta;
vcov(liwei.selection$linear.model);
X.new = rbind(c(1, -0.12), c(1, 0.12), c(1, 0.25));
y.hat.new = X.new %*% beta.WLS;
s2.y.hat.new = diag(X.new %*% vcov.beta %*% t(X.new));
s2.y.new = s2.y.hat.new + s2;
t.crit = qt(0.1, df = n - 2, lower.tail = F);
y.CI = cbind(y.hat.new - t.crit * sqrt(s2.y.hat.new), y.hat.new + t.crit * sqrt(s2.y.hat.new));
y.PI = cbind(y.hat.new - t.crit * sqrt(s2.y.new), y.hat.new + t.crit * sqrt(s2.y.new));
# CI and PI are consistent with predict.lm.
y.CI;
predict(liwei.selection$linear.model, newdata = list(x.bar = c(-0.12, 0.12, 0.25)),
    weights = c(1, 1, 1), interval = "confidence", level = 0.8);
y.PI;
predict(liwei.selection$linear.model, newdata = list(x.bar = c(-0.12, 0.12, 0.25)),
    weights = c(1, 1, 1), interval = "prediction", level = 0.8);



sunwei = list(
    average = c(16.62, 15.56, 4.45, 18.05, 39.09, 100.35, 83.66, -16.4, 20.95, 3.59,
        239.11, 63.93),
    fund = c(0.69, 7.16, 3.54, 3.69, 23.87, 136.11, 155.86, -11.36, 7.1, 0.87, 540.83,
        139.85),
    start = c(2020.84, 2020.67, 2020.6, 2020.45, 2020.2, 2018.84, 2018.55, 2017.9, 2016.75,
        2016.5, 2014.5, 2014.5),
    end = c(2021.8, 2021.8, 2021.6, 2021.8, 2021.8, 2021.54, 2021.8, 2019, 2020.21, 2018.73,
        2021.8, 2018.73)
);
length(sunwei$fund);
Plot.Manager(sunwei);
sunwei.selection = WLS.Select(sunwei);
summary(sunwei.selection$linear.model);
GLS.Plot.Family(sunwei, order = 1);



now = 2021.78;
taocan = list(
    name = "Tao Can",
    average = c(-1.64, 8.15, 48.34, 66.35, 67.55, 99.37, -6.8, 143.56, 287.32, 41.62),
    fund = c(-0.6, 11.94, 172.78, 111.67, 183.61, 166.42, 6.44, 99.5, 537.8, 36.93),
    start = c(2021.54, 2021.17, 2020.38, 2019.95, 2017.81, 2017.22, 2017.11, 2015.98,
              2014.27, 2011.52),
    end = c(now, now, now, now, now, now, 2019.06, 2021.74, now, 2016.31)
);
length(taocan$fund);
Plot.Manager(taocan);
taocan.selection = WLS.Select(taocan);
summary(taocan.selection$quadratic);
GLS.Plot.Family(taocan, order = 2);
WLS.Predict(taocan, order = 2);



now = 2021.78;
chenhao = list(
    average = c(-1.3, 3.25, 15.03, 18.34, 60.43, 45.86, 122.57, -1.4, 14.49, 154.78,
        33.72, 259.94, 255.2),
    fund = c(0.19, 11.83, 12.64, 25.59, 31.58, 88.16, 158.29, -18.79, -11.7, 318.2,
        33.57, 550.84, 392.14),
    start = c(2021.55, 2021.02, 2020.72, 2020.72, 2020.32, 2019.12, 2018.93, 2016.89,
        2015.43, 2015.04, 2014.88, 2014.36, 2012.73),
    end = c(now, now, 2021.02, now, now, 2020.58, now, 2018.97, 2020.43, now,
        2017.36, now, now)
);
chenhao.selection = WLS.Select(chenhao);
GLS.Plot.Family(chenhao, order = 2, moderated.only = T);



now = 2021.81;
xiaonan = list(
    average = c(12.31, 73.04, 95.54, 66.68, -0.07, 20.52, 66.44, -21.02, -1.87, 303.34),
    fund = c(-6.21, 22.67, 99.64, 143.2, 11, 41.7, 139.4, 40.6, 10.5, 452.44),
    start = c(2020.83, 2020.26, 2018.68, 2017.85, 2017.83, 2017.80, 2017.63, 2015.39,
        2015.21, 2012.74),
    end = c(now, now, now, now, 2019.23, 2020.43, 2021.41, 2016.6, 2016.6, now)
);
Plot.Manager(xiaonan);
xiaonan.selection = WLS.Select(xiaonan);
GLS.Plot.Family(xiaonan, order = 0);



now = 2021.81;
chenzhou = list(
    average = c(-0.65, 2.24, 41.28, 46.42, 59.41, 45.14, 69.33, 67.99, 124.95, 46.92,
        155, 249.56, 67.52, -16.16, -15.16, 4.44),
    fund = c(-0.3, 5.98, 35.93, 12.31, 95.59, 31.96, 32.37, 90, 138.9, 125.05, 177.94,
        448.6, 18.25, -13.2, -7.92, -21.08),
    start = c(2021.63, 2021, 2020.33, 2019.98, 2019.78, 2019.42, 2018.19, 2017.82, 2017.24,
        2017.21, 2015.07, 2013.78, 2012.98, 2010.04, 2009.96, 2008.26),
    end = c(now, now, now, 2021.53, now, 2020.57, now, now, now, 2020.57, now, now, 2017.07,
        2011.98, 2012.98, 2014.21)
);
chenzhou.selection = WLS.Select(chenzhou);
GLS.Plot.Family(chenzhou, order = 1);



now = 2021.82;
youlinfeng = list(
    name = "You Linfeng",
    average = c(8.92, 14.95, 83.2, 82.54, 8.02, 38.85, 3.88, 253.22, 290.41, 123.72, 91.49),
    fund = c(4.11, 17.59, 42.2, 63.2, 6.6, 42.4, 19.6, 36.8, 424.86, 222, 224.82),
    start = c(2020.93, 2018.82, 2017.23, 2017.19, 2016.66, 2016.05, 2015.93, 2014.44, 2012.21,
        2010.28, 2009.99),
    end = c(now, 2019.98, now, now, 2018.15, 2017.78, 2017.98, now, now, now, now)
);
youlinfeng.selection = WLS.Select(youlinfeng);
GLS.Plot.Family(youlinfeng, order = 0);



now = 2021.82;
wangzonghe = list(
    name = "Wang Zonghe",
    average = c(4.62, 14.53, 19.5, 21.09, 42.95, 46.97, 44.83, 80.6, 33.01, 71.33, 46.63,
        98.56, 15.97, 7.44, 18.42, 211.99, 74.07, 96.62, 213.86),
    fund = c(-23.63, 7.84, -18.56, -1.33, 15.33, 22.78, 47.36, 80.61, 5.87, 121.46, 110.95,
        132.18, 15.04, 6.9, 15.7, 287.2, 57.67, 55.52, 360.8),
    start = c(2021.03, 2020.73, 2020.73, 2020.42, 2020.3, 2020.23, 2020.11, 2019.32,
        2018.78, 2018.23, 2017.59, 2017.58, 2017.03, 2016.35, 2016.21, 2014.83, 2014.57,
        2012.33, 2010.91),
    end = c(now, 2021.03, now, now, now, now, now, now, 2020.16, now, 2020.61, now,
        2019.72, 2019.47, 2019.72, now, 2019.26, 2019.72, now)
);
wangzonghe.selection = WLS.Select(wangzonghe);
GLS.Plot.Family(wangzonghe, order = 1, moderated.only = T);

