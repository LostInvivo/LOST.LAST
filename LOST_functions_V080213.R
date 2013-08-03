

############################################################################################################################
getcontrol <- function(treatments)
{
	if(!require(tcltk)) stop('need tcltk')
    
	ntx <- length(treatments)
	tx.choice1 <- ''
    
   	tt<-tktoplevel()
    	tl<-tklistbox(tt,height=ntx,selectmode="single",background="white")
    	tkgrid(tklabel(tt,text="Select the Control Group" ))
    	tkgrid(tl)
    
    	for(i in (1:ntx))
    	{
      		tkinsert(tl,"end",treatments[i])
    	}
    
	OnOK<-function()
	{
		tx.choice1<<-treatments[as.integer(tkcurselection(tl))+1]
		tkdestroy(tt)
	}

	OK.but<-tkbutton(tt,text="    OK    ", command=OnOK)
	tkgrid(OK.but)
    
	tkwait.window(tt)
	return(tx.choice1)
}


############################################################################################################################
SplitVectorBy <- function(v,s,j){unlist(lapply( strsplit (v, s), function(x){x[[j]]}))} 
se <- function(x){sd(x)/sqrt(length(x))}
ci <- function(x){1.96*sd(x)/sqrt(length(x))}
cls <- function(){tmp<-flush.console()}

############################################################################################################################
Short2Long <- function(fname){
  tt <- read.csv(file=fname, head=FALSE, check.names =FALSE, skip=1, stringsAsFactors=FALSE)
  h <- read.csv(file=fname, head=FALSE, nrows=1, stringsAsFactors=FALSE)
  h[1,1] <- "time"
  colnames(tt) <- c(h[1,])
  cmpds <- colnames(tt)[-1]
  tabc <- table(cmpds)
  #remove time
  ttc <- colnames(tt)[-1]
  time <- tt[,1]
  tt <- tt[,-1]
  
  tW <- data.frame("1","1",1,"1");
  colnames(tW) <- c("time", "aid","value","treatid")
  tW <- tW[-1,]
  
  for (tr in names(tabc)){
    colids <- which(ttc %in% tr)
    t0 <- tt[,colids] 
    colnames(t0) <- paste("an",colids,sep="")
    t0$time <- time
    t1 <- melt(t0, id=c("time"))
    t1$treatid <- tr
    colnames(t1) <- c("time", "aid","value","treatid")
    tW <- rbind(tW, t1)
  }
  tW
}

############################################################################################################################
plot_time_dose <- function( tR, title="", xlab="", ylab=""){
  xr <- range(tR$time)
  df.plt <- summaryBy( value ~ time + treatid, tR, FUN=c(mean,se))
  df.plt$time<-as.numeric(df.plt$time)
  ggplot(df.plt, aes(x=time, y=value.mean, colour=treatid, group=treatid)) + 
    geom_errorbar(aes(ymin=value.mean-value.se, ymax=value.mean+value.se), width=5,size=1.) + 
    geom_line(size=1.5) + geom_point(aes(symbol=treatid)) + theme_bw() + scale_x_continuous(breaks=seq(xr[1],xr[2],10)) +
    theme(text = element_text(size=18)) + xlab(xlab) + ylab(ylab)+ ggtitle(title) 
}

############################################################################################################################
# this is fix of lsmeans to make it work with newer version of nlme
# patched by Bob Gerwien 07.11.13

lsmeans2<-
function (object, specs, adjust = c("auto", "tukey", "sidak", 
    "scheffe", p.adjust.methods), conf = 0.95, at, trend, contr = list(), 
    cov.reduce = function(x, name) mean(x), fac.reduce = function(coefs, 
        lev) apply(coefs, 2, mean), glhargs = NULL, lf = FALSE, 
    ...) 
{
    if (missing(specs)) 
        stop("Must specify specs, e.g. 'pairwise ~ treatment'")
    if (!is.null(glhargs)) {
        if (!require("multcomp")) {
            glhargs = NULL
            warning("'glhargs' option disabled because 'multcomp' package not installed")
        }
        else {
            if (!is.null(glhargs$df)) 
                glhargs$df = as.integer(max(1, 0.2 + glhargs$df))
        }
    }
    if (is.logical(cov.reduce)) {
        if (cov.reduce) 
            cov.reduce = function(x, name) mean(x)
        else cov.reduce = function(x, name) sort(unique(x))
    }
    adjtbl = c("auto", "tukey", "sidak", "scheffe", p.adjust.methods)
    no.adj = pmatch("none", adjtbl)
    adj = pmatch(adjust, adjtbl)[1]
    if (is.na(adj)) {
        adj = 1
        warning("Unknown or non-unique `adjust' method -- automatic method will be used")
    }
    autoadj = (adj == 1)
    trend.flag = !missing(trend)
#    if (inherits(object, "gls")) 
#        Terms = getCovariateFormula(object)
#    else Terms = terms(object)
    Terms = terms(object)
    formrhs = formula(Terms)
    ddfm = adjV = NULL
    if (inherits(object, "mer") || inherits(object, "merMod")) {
        if (!isLMM(object) && !isGLMM(object)) 
            stop("Can't handle a nonlinear mixed model")
        thecall = slot(object, "call")
        bhat = fixef(object)
        contrasts = attr(model.matrix(object), "contrasts")
        if (isLMM(object)) {
            if (require("pbkrtest")) {
                adjV = vcovAdj(object, 0)
                ddfm = function(k, se) .KRdf.mer(adjV, V, k, 
                  se * se)
            }
            else warning("Install package 'pbkrtest' to obtain bias corrections and degrees of freedom")
        }
    }
    else if (inherits(object, "lme")) {
        thecall = object$call
        bhat = fixef(object)
        contrasts = object$contrasts
    }
    else if (inherits(object, "gls")) {
        thecall = object$call
        bhat = coef(object)
        contrasts = object$contrasts
        the.df = object$dims$N - object$dims$p
        ddfm = function(k, se) the.df
    }
    else if (inherits(object, "lm")) {
        thecall = object$call
        bhat = coef(object)
        contrasts = attr(model.matrix(object), "contrasts")
        if (!(family(object)$family %in% c("binomial", "poisson"))) 
            if (!is.na(object$df.residual)) 
                ddfm = function(k, se) object$df.residual
    }
    else stop(paste("Can't handle an object of class", class(object)[1]))
    if (is.null(adjV)) 
        V = vcov(object)
    else V = adjV
    used = which(!is.na(bhat))
    not.used = which(is.na(bhat))
    bhat = bhat[used]
    if (length(bhat) != nrow(V)) 
        stop("Something's wrong -- Mismatch between vcov() and non-missing coef() results")
    null.basis = NULL
    if (length(not.used) > 0) {
        tR = t(qr.R(object$qr))
        if (ncol(tR) < nrow(tR)) 
            tR = cbind(tR, matrix(0, nrow = nrow(tR), ncol = nrow(tR) - 
                ncol(tR)))
        rank = object$qr$rank
        for (i in (rank + 1):nrow(tR)) tR[i, i] = 1
        null.basis = qr.resid(qr(tR[, seq_len(rank)]), tR[, -seq_len(rank)])
        if (!is.matrix(null.basis)) 
            null.basis = matrix(null.basis, ncol = 1)
        null.basis[object$qr$pivot, ] = null.basis
    }
    nm = all.vars(formrhs)
    anm = all.names(formrhs)
    coerced = anm[1 + grep("factor|ordered", anm)]
    form = as.formula(paste("~", paste(nm, collapse = "+")))
    envir = attr(Terms, ".Environment")
    X = model.frame(form, eval(thecall$data, envir = envir), 
        subset = eval(thecall$subset, enclos = envir), na.action = na.omit, 
        drop.unused.levels = TRUE)
    baselevs = xlev = matdat = list()
    if (is.character(specs)) 
        specs = as.list(specs)
    if (!is.list(specs)) 
        specs = list(specs)
    all.var.names = names(X)
    for (xname in all.var.names) {
        obj = X[[xname]]
        if (is.factor(obj)) {
            xlev[[xname]] = levels(obj)
            if (!missing(at) && !is.null(at[[xname]])) 
                baselevs[[xname]] = at[[xname]]
            else baselevs[[xname]] = levels(obj)
        }
        else if (is.matrix(obj)) {
            matdat[[xname]] = apply(obj, 2, cov.reduce, xname)
        }
        else {
            if (length(grep(xname, coerced)) > 0) 
                baselevs[[xname]] = sort(unique(obj))
            else {
                if (!missing(at) && !is.null(at[[xname]])) 
                  baselevs[[xname]] = at[[xname]]
                else baselevs[[xname]] = cov.reduce(obj, xname)
            }
        }
    }
    if (trend.flag) {
        if (!is.character(trend)) 
            stop("'trend' must be of character type")
        trend.xnm = trend[1]
        trend.x = X[[trend.xnm]]
        if (is.null(trend.x)) {
            trend.xnm = try(all.vars(as.formula(paste("~", trend))), 
                silent = TRUE)
            if (inherits(trend.xnm, "try-error")) 
                trend.xnm = character(0)
            trend.h = -1
        }
        else {
            if (!is.numeric(trend.x)) 
                stop("'trend' must refer to a numeric predictor")
            trend.h = diff(range(trend.x)) * 0.001
            baselevs[[trend]] = baselevs[[trend]][1] + c(-1, 
                1) * trend.h/2
            sidx = match(trend, names(baselevs))
            baselevs = c(baselevs[sidx], baselevs[-sidx])
        }
    }
    yidx = attr(Terms, "response")
    if (yidx > 0) {
        yname = as.character(attr(Terms, "variables")[[1 + yidx]])
        if (!is.na(match(yname, names(baselevs))[1])) 
            baselevs[[yname]] = NA
    }
    grid = do.call(expand.grid, baselevs)
    for (nm in names(matdat)) grid[[nm]] = matrix(rep(matdat[[nm]], 
        each = nrow(grid)), nrow = nrow(grid))
    m = model.frame(Terms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(Terms, m, contrasts.arg = contrasts)
    if (trend.flag) {
        if (trend.h < 0) {
            term.nm = gsub(" ", "", dimnames(X)[[2]])
            trend = gsub(" ", "", trend)
            trend.idx = match(trend, term.nm)
            if (is.na(trend.idx)) 
                stop(paste("trend value '", trend, "' is neither a variable nor a model term", 
                  sep = ""))
            prev.X = X
            for (i in seq_len(ncol(X))) {
                trm = term.nm[i]
                trm.pieces = strsplit(trm, ":")[[1]]
                trm.mat = match(trend, trm.pieces)
                if (is.na(trm.mat)) 
                  X[, i] = 0
                else {
                  trm.otr = paste(trm.pieces[-trm.mat], collapse = ":")
                  trm.ref = match(trm.otr, term.nm)
                  if (is.na(trm.ref)) 
                    X[, i] = 0 + (i == trend.idx)
                  else X[, i] = prev.X[, trm.ref]
                }
            }
        }
        else {
            evens = 2 * (seq_len(nrow(X)/2))
            X = (X[evens, ] - X[evens - 1, ])/trend.h
            baselevs[[trend]] = baselevs[[trend]][1] + trend.h/2
            grid = grid[evens, , drop = FALSE]
        }
    }
    if (length(coerced) > 0) 
        grid = do.call("expand.grid", baselevs)
    allFacs = all.var.names
    row.indexes = array(seq_len(nrow(X)), sapply(baselevs, length))
    mod.terms = attr(Terms, "term.labels")
    some.term.contains = function(facs) {
        if (trend.flag) 
            facs = union(facs, trend.xnm)
        for (trm in mod.terms) {
            flag = all(sapply(facs, function(f) length(grep(f, 
                trm)) > 0))
            if (flag) 
                if (length(all.vars(as.formula(paste("~", trm)))) > 
                  length(facs)) 
                  return(TRUE)
        }
        return(FALSE)
    }
    results = list()
    for (i in seq_len(length(specs))) {
        form = specs[[i]]
        if (is.character(form)) 
            form = as.formula(paste("~", form))
        if (!inherits(form, "formula")) 
            stop(paste("Incorrect formula specification:", form))
        method = byfacs = NULL
        if (length(form) == 3) {
            method = all.vars(form[[2]])[1]
            form = form[-2]
        }
        facs = all.vars(form)
        facs.lbl = paste(facs, collapse = ":")
        if (some.term.contains(facs)) 
            warning(paste("lsmeans of", facs.lbl, "may be misleading due to interaction with other predictor(s)"))
        ln = if (any(sapply(facs, function(nm) length(grep(nm, 
            allFacs)) == 0))) 
            stop(paste("Unknown factor(s) in specification:", 
                paste(form, collapse = " ")))
        b = strsplit(as.character(form[2]), "\\|")[[1]]
        if (length(b) > 1) 
            byfacs = all.vars(as.formula(paste("~", b[2])))
        levs = list()
        for (f in facs) levs[[f]] = baselevs[[f]]
        combs = do.call("expand.grid", levs)
        RI = plyr:::splitter_a(row.indexes, match(facs, names(baselevs)))
        K = sapply(RI, function(idx) {
            fac.reduce(X[idx, , drop = FALSE], "")
        })
        rnames = dimnames(K)[[2]] = apply(combs, 1, paste, collapse = ", ")
        do.est = function(k) {
            est = se = df = NA
            estimable = TRUE
            if (!is.null(null.basis)) {
                estimable = all(abs(apply(null.basis, 2, function(x) sum(k * 
                  x))) < 1e-04)
            }
            if (estimable) {
                k = k[used]
                est = sum(k * bhat)
                se = sqrt(sum(k * (V %*% k)))
                if (!is.null(ddfm)) 
                  df = ddfm(k, se)
            }
            c(estimate = est, SE = se, df = df)
        }
        adj.p.value = function(t, df, meth, fam.size, n.contr) {
            abst = abs(t)
            if (meth <= 4) 
                switch(meth, NA, ptukey(sqrt(2) * abst, fam.size, 
                  zapsmall(df), lower.tail = FALSE), 1 - (1 - 
                  2 * pt(abst, df, lower.tail = FALSE))^n.contr, 
                  pf(t^2/(fam.size - 1), fam.size - 1, df, lower.tail = FALSE))
            else p.adjust(2 * pt(abst, df, lower.tail = FALSE), 
                adjtbl[meth], n = n.contr)
        }
        if (!trend.flag) {
            effname = "lsmean"
            lsmentry = paste(facs.lbl, "lsmeans")
        }
        else {
            effname = paste(trend, "trend", sep = ".")
            lsmentry = paste(effname, "by", facs.lbl)
        }
        if (lf) {
            results[[lsmentry]] = t(K)
        }
        else {
            lsms = as.data.frame(t(apply(K, 2, do.est)))
            names(lsms)[1] = effname
            lsms = cbind(combs, lsms)
            if (conf > 1) 
                conf = conf/100
            if ((conf < 1) && (conf > 0.01)) {
                if (is.null(ddfm)) {
                  me = qnorm((1 - conf)/2, lower.tail = FALSE) * 
                    lsms$SE
                  lsms$asymp.LCL = lsms[[effname]] - me
                  lsms$asymp.UCL = lsms[[effname]] + me
                }
                else {
                  me = qt((1 - conf)/2, lsms$df, lower.tail = FALSE) * 
                    lsms$SE
                  lsms$lower.CL = lsms[[effname]] - me
                  lsms$upper.CL = lsms[[effname]] + me
                }
            }
            attr(lsms, "print.row.names") = FALSE
            class(lsms) = c("data.frame.lsm", "data.frame")
            results[[lsmentry]] = lsms
        }
        if (!is.null(method)) {
            fn = paste(method, "lsmc", sep = ".")
            confcn = if (exists(fn, mode = "function")) 
                get(fn)
            else NULL
            if (is.null(byfacs)) 
                bylist = list(seq_len(nrow(combs)))
            else {
                bg = list()
                for (f in byfacs) bg[[f]] = baselevs[[f]]
                bygrid = do.call("expand.grid", bg)
                bylist = lapply(seq_len(nrow(bygrid)), function(row) {
                  bylevs = bygrid[row, ]
                  if (length(byfacs) > 1) 
                    flags = apply(combs[, byfacs], 1, function(r) all(r == 
                      bylevs))
                  else flags = combs[, byfacs] == bylevs
                  which(flags)
                })
                bylabs = apply(bygrid, 1, paste, collapse = ",")
                bycols = sapply(byfacs, grep, names(combs))
                rnames = combs[, -bycols]
                if (!is.null(ncol(rnames))) 
                  rnames = apply(rnames, 1, paste, collapse = ",")
            }
            Clist = list()
            zer = rep(0, ncol(K))
            nby = length(bylist)
            for (i in seq_len(nby)) {
                rows = bylist[[i]]
                cl = if (is.null(confcn)) 
                  contr[[method]]
                else confcn(rnames[rows], ...)
                if (is.null(cl)) 
                  stop(paste("Unknown contrast family:", method))
                clx = lapply(cl, function(cc) {
                  if (length(cc) != length(rows)) 
                    stop(paste(length(cc), " contrast coefficients in '", 
                      method, "' when ", length(rows), " were expected", 
                      sep = ""))
                  ccc = zer
                  ccc[rows] = cc
                  ccc
                })
                if (nby > 1) 
                  names(clx) = paste(names(clx), "|", bylabs[i])
                Clist = c(Clist, clx)
            }
            methdesc = attr(cl, "desc")
            if (!is.null(null.basis) && !is.null(glhargs)) {
                warning("Error may occur in 'glht' due to rank deficiency")
            }
            if (lf || !is.null(glhargs)) {
                KK = t(sapply(Clist, function(con) {
                  nz = which(abs(con) > 1e-04)
                  K[, nz] %*% con[nz]
                }))
                if (lf) {
                  dimnames(KK)[[2]] = row.names(K)
                  ctbl = KK
                }
                else {
                  args = c(list(model = object, linfct = KK[, 
                    used]), glhargs)
                  ctbl = summary(do.call("glht", args))
                }
            }
            else {
                if (is.null(methdesc)) 
                  methdesc = method
                adjattr = attr(cl, "adjust")
                if (autoadj) 
                  adj = ifelse(is.null(adjattr), no.adj, pmatch(adjattr, 
                    adjtbl))
                if (is.na(adj)) 
                  adj = no.adj
                ctbl = as.data.frame(t(sapply(Clist, function(con) {
                  nz = which(abs(con) > 1e-04)
                  k = K[, nz] %*% con[nz]
                  do.est(k)
                })))
                n.fam = nrow(lsms)/nby
                n.contr = sum(!is.na(ctbl$estimate))
                if (!is.null(ddfm)) {
                  ctbl$t.ratio = round(ctbl$estimate/ctbl$SE, 
                    5)
                  ctbl$p.value = round(adj.p.value(ctbl$t.ratio, 
                    ctbl$df, adj, n.fam, n.contr), 5)
                }
                else {
                  ctbl$z.ratio = round(ctbl$estimate/ctbl$SE, 
                    5)
                  ctbl$p.value = round(adj.p.value(ctbl$z.ratio, 
                    10000, adj, n.fam, n.contr), 5)
                }
                attr(ctbl, "mesg") = if (adj == 2) 
                  paste("p values are adjusted using the", adjtbl[adj], 
                    "method for", n.fam, "means")
                else if (adj < length(adjtbl)) 
                  paste("p values are adjusted using the", adjtbl[adj], 
                    "method for", n.contr, "tests")
                else "p values are not adjusted"
                attr(ctbl, "print.row.names") = TRUE
                class(ctbl) = c("data.frame.lsm", "data.frame")
            }
            results[[paste(facs.lbl, methdesc)]] = ctbl
        }
    }
    if (!lf) 
        class(results) = c("lsm", "list")
    results
}

############################################################################################################################



