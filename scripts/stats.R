scale_01 = function(x) {
  (x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
}

scale_01_mat = function(m) {
  (m - min(c(m), na.rm=T)) / (max(c(m)) - min(c(m)))
}

calc_geom_mean = function(x) {
  x = na.omit(x)
  prod(x) ^ (1 / length(x))
}

calc_cv = function(x) {
  sd(x) / abs(mean(x))
}

calc_mean_sd = function(x, name) {
  d = data.frame(mean(x,na.rm=T), sd(x,na.rm=T))
  colnames(d) = paste(name, c("mean", "sd"), sep="_")
  d
}

calc_kl = function(x,y) {
  sum(x * log2(x/y))
}

calc_jsd = function(x,y) {
  x = x[x>0]
  y = y[y>0]
  m = (x + y) / 2
  sqrt(.5 * calc_kl(x,m) + .5 * calc_kl(y,m))
}

calc_eu = function(x, y) {
  sqrt(sum((x - y)^2))
}

sem = function(x) {
  x = na.omit(x)
  sd(x)/sqrt(length(x))
}

mean_na = function(x) {
  x = x[is.finite(x)]
  mean(x, na.rm=T)
}

median_na = function(x) {
  x = x[is.finite(x)]
  median(x, na.rm=T)
}

iqr = function(x) {
  qs = quantile(na.omit(x), probs = c(.25, .75))
  diff(qs)
}

sd_na = function(x) {
  x = x[is.finite(x)]
  sd(x, na.rm=T)
}

sem_na = function(x) {
  x = na.omit(x)
  sd(x)/length(x)
}

#' https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
#' bootstrap validation 
weighted_sem = function(x, w) {
  na_ind = is.na(x)
  x = x[!na_ind]
  w = w[!na_ind]
  w1 = (w - min(w)) / (max(w) - min(w))
  ws = sum(w1)
  n = length(x)
  std_var = (n / ((n-1) * ws^2)) * sum(w1^2 * (x - mean(x))^2)
  sqrt(std_var)
}

var_na = function(x) {
  x = x[is.finite(x)]
  var(x, na.rm=T)
}
q75 = function(x) {
  x = na.omit(x)
  quantile(x, probs=c(.75))
}

zscore_clean = function(d) {
  d = na.omit(d)
  (d - mean(d)) / sd(d)
}
zscore_my = function(d, drel=NULL, default_sd=1, baseline=NULL) {
  if (is.null(drel) & !is.null(baseline)) {
    drel = d[baseline]
  }
  if (length(drel) < 2) {
    sd_value = default_sd
  } else {
    sd_value = sd(drel, na.rm=T)    
  }
  mean_value = mean(drel, na.rm=T)
  (d - mean_value) / sd_value
}

frac_to_perc = function(x) {
  x * 100
}

perc_change = function(d, drel=NULL, baseline=NULL) {
  if (is.null(drel) & !is.null(baseline)) {
    drel = d[baseline]
  }
  baseline_mean =  mean(drel, na.rm=T)
  100 * (d - baseline_mean) / baseline_mean
}

calc_diff = function(d, drel=NULL, baseline=NULL) {
  if (is.null(drel) & !is.null(baseline)) {
    drel = d[baseline]
  }
  baseline_mean =  mean(drel, na.rm=T)
  d - baseline_mean
}

#' Tissue specificty score
#'
#' @param x 
#'
#' @return scalar
#' @export
#' @references http://www.sciencedirect.com/science/article/pii/S0092867410000796
#' @examples
tissue_specificity = function(x, scale=F) {
  if (scale)
    x = (x - (min(x))) / (max(x) - (min(x)))
  x = x[x>0]
  x_norm = x /sum(x)
  #sum(x_norm * log(x_norm / (1/length(x))))
  sum(x_norm * log(x_norm / mean(x_norm)))
}

calc_tissue_specificity= function(x, scale=F) {
  if (scale)
    x = (x - (min(x))) / (max(x) - (min(x)))
  x_norm = x /sum(x)
  y = x_norm * log(x_norm / (1/length(x)))
  y[is.nan(y)] = 0
  sum(y)
}

calc_gene_specificity = function(x, scale=F) {
  if (scale)
    x = (x - (min(x))) / (max(x) - (min(x)))
  #x = x[x>0]
  x_norm = x /sum(x)
  x_norm * log(x_norm / mean(x_norm) + .01)
}

# calc_tau = function(x) {
#   x_hat = x / max(x)
#   tau = 
# }



