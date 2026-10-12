library(ggplot2)

TAU      <- 2
N        <- 10000
ISMETLES <- 5000
SEED     <- 2026

# kimenet a szkript melletti eredmenyek mappaba, Rscript es RStudio alatt is
szkript_mappa <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("^--file=", "", arg[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(p))
  }
  getwd()
}
KIMENET <- file.path(szkript_mappa(), "eredmenyek")


# adatgeneralo folyamatok

normalis_df <- function(n, nevek) {
  df <- as.data.frame(matrix(rnorm(n * length(nevek)), nrow = n))
  names(df) <- nevek
  df
}

kezeles_es_kimenet <- function(df, index, y0) {
  df$e_valodi <- plogis(index)
  df$kezeles  <- rbinom(nrow(df), 1, df$e_valodi)
  df$Y        <- y0 + TAU * df$kezeles + rnorm(nrow(df))
  df
}

# a: kezeles logit-linearis a megfigyelt confounderekben
dgp_mukodik <- function(n, szelekcio = 1) {
  df    <- normalis_df(n, c("X1", "X2", "X3"))
  index <- szelekcio * (0.8 * df$X1 - 0.6 * df$X2 + 0.5 * df$X3)
  y0    <- 1.5 * df$X1 - 1.0 * df$X2 + 1.0 * df$X3
  kezeles_es_kimenet(df, index, y0)
}

# b: kezeles es kimenet negyzetesen fugg X1, X2-tol, a linearis logit ezt nem kapja el
dgp_rossz_specifikacio <- function(n) {
  df <- normalis_df(n, c("X1", "X2", "X3"))
  df$X1negyzet <- df$X1^2
  df$X2negyzet <- df$X2^2
  index <- 0.6 * (df$X1negyzet - 1) + 0.6 * (df$X2negyzet - 1) + 0.3 * df$X3
  y0    <- 1.2 * df$X1negyzet + 1.2 * df$X2negyzet + 1.0 * df$X3
  kezeles_es_kimenet(df, index, y0)
}

# c: U hat a kezelesre es a kimenetre is, de nem figyeljuk meg
dgp_rejtett_confounder <- function(n) {
  df    <- normalis_df(n, c("X1", "X2", "U"))
  index <- 0.6 * df$X1 - 0.5 * df$X2 + 1.0 * df$U
  y0    <- 1.0 * df$X1 - 1.0 * df$X2 + 1.5 * df$U
  kezeles_es_kimenet(df, index, y0)
}

# d: mint a, csak negyszer erosebb szelekcio
dgp_nincs_overlap <- function(n) dgp_mukodik(n, szelekcio = 4)

# e: Z erosen hat a kezelesre, a kimenetre nem
dgp_instrumentum <- function(n) {
  df    <- normalis_df(n, c("X1", "X2", "Z"))
  index <- 0.6 * df$X1 - 0.5 * df$X2 + 3.0 * df$Z
  y0    <- 1.0 * df$X1 - 1.0 * df$X2
  kezeles_es_kimenet(df, index, y0)
}


# becslok

becsult_ps <- function(df, valtozok) {
  modell <- glm(reformulate(valtozok, response = "kezeles"), family = binomial, data = df)
  unname(fitted(modell))
}

# minden honnan elemhez a legkozelebbi hova elem indexe
legkozelebbi <- function(honnan, hova) {
  rend <- order(hova)
  s    <- hova[rend]
  pos  <- pmin(pmax(findInterval(honnan, s), 1), length(s) - 1)
  balra_kozelebb <- (honnan - s[pos]) <= (s[pos + 1] - honnan)
  rend[ifelse(balra_kozelebb, pos, pos + 1)]
}

# 1:1 legkozelebbi szomszed a logit ps-en, visszatevessel, mindket iranyban (ATE)
psm <- function(df, ps, caliper = NULL) {
  lp       <- qlogis(pmin(pmax(ps, 1e-6), 1 - 1e-6))
  kezelt   <- which(df$kezeles == 1)
  kontroll <- which(df$kezeles == 0)

  par <- integer(nrow(df))
  par[kezelt]   <- kontroll[legkozelebbi(lp[kezelt], lp[kontroll])]
  par[kontroll] <- kezelt[legkozelebbi(lp[kontroll], lp[kezelt])]

  egyeni_hatas <- ifelse(df$kezeles == 1, df$Y - df$Y[par], df$Y[par] - df$Y)
  megtart <- rep(TRUE, nrow(df))
  if (!is.null(caliper)) megtart <- abs(lp - lp[par]) <= caliper * sd(lp)
  list(ate = mean(egyeni_hatas[megtart]), par = par, megtart = megtart)
}

naiv_kulonbseg <- function(df) {
  mean(df$Y[df$kezeles == 1]) - mean(df$Y[df$kezeles == 0])
}

auc <- function(kezeles, ps) {
  n1 <- sum(kezeles == 1)
  n0 <- sum(kezeles == 0)
  (sum(rank(ps)[kezeles == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# standardizalt atlagkulonbseg, par nelkul a nyers, par-ral a parositott mintan
smd <- function(df, valtozok, par = NULL, megtart = NULL) {
  X <- as.matrix(df[valtozok])
  k <- df$kezeles == 1
  szoras <- sqrt((apply(X[k, , drop = FALSE], 2, var) + apply(X[!k, , drop = FALSE], 2, var)) / 2)
  if (is.null(par)) {
    kulonbseg <- colMeans(X[k, , drop = FALSE]) - colMeans(X[!k, , drop = FALSE])
  } else {
    i <- seq_len(nrow(df))
    kezelt_oldal   <- ifelse(k, i, par)[megtart]
    kontroll_oldal <- ifelse(k, par, i)[megtart]
    kulonbseg <- colMeans(X[kezelt_oldal, , drop = FALSE]) - colMeans(X[kontroll_oldal, , drop = FALSE])
  }
  kulonbseg / szoras
}


# forgatokonyvek, az elso psm valtozat a vizsgalt, a masodik a referencia

FORGATOKONYVEK <- list(
  list(
    kod = "A", nev = "jól specifikált PS-modell", dgp = dgp_mukodik,
    psm = list(list(cimke = "PSM", valtozok = c("X1", "X2", "X3"))),
    balance = c(X1 = "X1", X2 = "X2", X3 = "X3")
  ),
  list(
    kod = "B", nev = "rosszul specifikált PS-modell", dgp = dgp_rossz_specifikacio,
    psm = list(
      list(cimke = "PSM, lineáris logit", valtozok = c("X1", "X2", "X3")),
      list(cimke = "PSM, helyes specifikáció", valtozok = c("X1", "X2", "X3", "X1negyzet", "X2negyzet"))
    ),
    balance = c(X1 = "X1", X2 = "X2", X3 = "X3", X1negyzet = "X1² (kimaradt)", X2negyzet = "X2² (kimaradt)")
  ),
  list(
    kod = "C", nev = "meg nem figyelt confounder", dgp = dgp_rejtett_confounder,
    psm = list(
      list(cimke = "PSM, U nélkül", valtozok = c("X1", "X2")),
      list(cimke = "PSM, U-val", valtozok = c("X1", "X2", "U"))
    ),
    balance = c(X1 = "X1", X2 = "X2", U = "U (rejtett)")
  ),
  list(
    kod = "D", nev = "nincs overlap", dgp = dgp_nincs_overlap,
    psm = list(
      list(cimke = "PSM", valtozok = c("X1", "X2", "X3")),
      list(cimke = "PSM, caliper 0,2 szórás", valtozok = c("X1", "X2", "X3"), caliper = 0.2)
    ),
    balance = c(X1 = "X1", X2 = "X2", X3 = "X3")
  ),
  list(
    kod = "E", nev = "instrumentum a PS-modellben", dgp = dgp_instrumentum,
    psm = list(
      list(cimke = "PSM, Z-vel", valtozok = c("X1", "X2", "Z")),
      list(cimke = "PSM, Z nélkül", valtozok = c("X1", "X2"))
    ),
    balance = c(X1 = "X1", X2 = "X2", Z = "Z (instrumentum)")
  )
)


# monte carlo

monte_carlo <- function(ismetles = ISMETLES, n = N, seed = SEED) {
  set.seed(seed)
  eredmeny <- list()
  for (fk in FORGATOKONYVEK) {
    cat(sprintf("  %s  %s ...\n", fk$kod, fk$nev))
    cimkek  <- c("naiv különbség", sapply(fk$psm, function(v) v$cimke))
    becsles <- matrix(NA_real_, ismetles, length(cimkek))
    auc_m   <- becsles
    megtart <- becsles
    for (r in seq_len(ismetles)) {
      df <- fk$dgp(n)
      becsles[r, 1] <- naiv_kulonbseg(df)
      megtart[r, 1] <- 1
      for (j in seq_along(fk$psm)) {
        v  <- fk$psm[[j]]
        ps <- becsult_ps(df, v$valtozok)
        m  <- psm(df, ps, v$caliper)
        becsles[r, j + 1] <- m$ate
        auc_m[r, j + 1]   <- auc(df$kezeles, ps)
        megtart[r, j + 1] <- mean(m$megtart)
      }
    }
    eredmeny[[fk$kod]] <- data.frame(
      kod        = fk$kod,
      becslo     = rep(cimkek, each = ismetles),
      becsles    = as.vector(becsles),
      auc        = as.vector(auc_m),
      megtartott = as.vector(megtart)
    )
  }
  do.call(rbind, c(eredmeny, make.row.names = FALSE))
}

osszefoglalo <- function(mc) {
  kulcs  <- paste(mc$kod, mc$becslo, sep = "|")
  reszek <- split(mc, factor(kulcs, levels = unique(kulcs)))
  nevek  <- setNames(sapply(FORGATOKONYVEK, function(fk) fk$nev), sapply(FORGATOKONYVEK, function(fk) fk$kod))
  tabla <- do.call(rbind, lapply(reszek, function(d) data.frame(
    kod              = d$kod[1],
    forgatokonyv     = nevek[[d$kod[1]]],
    becslo           = d$becslo[1],
    atlag            = mean(d$becsles),
    torzitas         = mean(d$becsles) - TAU,
    szoras           = sd(d$becsles),
    rmse             = sqrt(mean((d$becsles - TAU)^2)),
    auc              = mean(d$auc),
    megtartott_arany = mean(d$megtartott)
  )))
  rownames(tabla) <- NULL
  tabla
}


# abrak

HATTER <- "#fcfcfb"; TINTA <- "#0b0b0b"; HALVANY <- "#52514e"; RACS <- "#e4e3df"
KEK <- "#2a78d6"; NARANCS <- "#eb6834"; SZURKE <- "#8a8983"

tema <- theme_minimal(base_size = 12) +
  theme(
    plot.background     = element_rect(fill = HATTER, colour = NA),
    plot.title          = element_text(face = "bold", size = 14, colour = TINTA, margin = margin(b = 10)),
    plot.title.position = "plot",
    plot.margin         = margin(12, 16, 10, 12),
    panel.grid.minor    = element_blank(),
    panel.grid.major    = element_line(colour = RACS, linewidth = 0.3),
    axis.title          = element_text(colour = HALVANY, size = 11),
    axis.title.x        = element_text(margin = margin(t = 8)),
    axis.title.y        = element_text(margin = margin(r = 8)),
    axis.text           = element_text(colour = HALVANY, size = 10),
    legend.position     = "bottom",
    legend.title        = element_blank(),
    legend.text         = element_text(size = 11),
    legend.key.spacing.x = unit(14, "pt")
  )

abracim <- function(fk) paste(fk$nev, "eset")

abra_ps_eloszlas <- function(fk, minta) {
  adat <- data.frame(
    ps      = minta$ps,
    csoport = factor(ifelse(minta$df$kezeles == 1, "kezelt", "kontroll"), levels = c("kontroll", "kezelt"))
  )
  hatarok <- seq(0, 1, by = 0.025)
  ggplot(adat, aes(x = ps)) +
    geom_histogram(aes(fill = csoport), breaks = hatarok, position = "identity", alpha = 0.45) +
    stat_bin(aes(colour = csoport), geom = "step", breaks = hatarok, position = "identity",
             direction = "mid", pad = TRUE, linewidth = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = c(kontroll = KEK, kezelt = NARANCS)) +
    scale_colour_manual(values = c(kontroll = KEK, kezelt = NARANCS)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = abracim(fk), x = "becsült propensity score", y = "megfigyelések száma") +
    tema + theme(panel.grid.major.x = element_blank())
}

abra_balance <- function(fk, minta) {
  valtozok <- names(fk$balance)
  adat <- data.frame(
    valtozo = factor(unname(fk$balance), levels = rev(unname(fk$balance))),
    elotte  = abs(smd(minta$df, valtozok)),
    utana   = abs(smd(minta$df, valtozok, minta$par, minta$megtart))
  )
  pontok <- rbind(
    data.frame(valtozo = adat$valtozo, smd = adat$elotte, mikor = "párosítás előtt"),
    data.frame(valtozo = adat$valtozo, smd = adat$utana,  mikor = "párosítás után")
  )
  ggplot(adat, aes(y = valtozo)) +
    geom_segment(aes(x = elotte, xend = utana, yend = valtozo), colour = RACS, linewidth = 0.8) +
    geom_vline(xintercept = 0.1, linetype = "dashed", colour = HALVANY, linewidth = 0.4) +
    geom_point(data = pontok, aes(x = smd, colour = mikor, fill = mikor), shape = 21, size = 3, stroke = 1) +
    scale_colour_manual(values = c("párosítás előtt" = SZURKE, "párosítás után" = KEK)) +
    scale_fill_manual(values = c("párosítás előtt" = HATTER, "párosítás után" = KEK)) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0.04, 0.06))) +
    labs(title = abracim(fk), x = "standardizált átlagkülönbség (abszolút érték)", y = NULL) +
    tema + theme(panel.grid.major.y = element_blank())
}

abra_becslesek <- function(fk, mc) {
  adat   <- mc[mc$kod == fk$kod, ]
  cimkek <- unique(adat$becslo)
  suruseg <- do.call(rbind, lapply(cimkek, function(cimke) {
    s <- density(adat$becsles[adat$becslo == cimke], n = 300)
    data.frame(becslo = cimke, x = s$x, y = s$y)
  }))
  suruseg$becslo <- factor(suruseg$becslo, levels = cimkek)
  szinek <- setNames(c(SZURKE, KEK, NARANCS)[seq_along(cimkek)], cimkek)

  ggplot(suruseg, aes(x = x, y = y)) +
    geom_area(aes(fill = becslo), position = "identity", alpha = 0.18) +
    geom_line(aes(colour = becslo), linewidth = 0.9) +
    geom_vline(xintercept = TAU, linetype = "dashed", colour = TINTA, linewidth = 0.5) +
    scale_colour_manual(values = szinek) +
    scale_fill_manual(values = szinek) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = abracim(fk), x = "ATE becslés", y = NULL) +
    tema + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
}

ment <- function(abra, fk, nev, magassag = 4.2) {
  mappa <- file.path(KIMENET, fk$kod)
  dir.create(mappa, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(mappa, sprintf("%s_%s.png", fk$kod, nev)), abra,
         width = 6.5, height = magassag, dpi = 200)
}


# futtatas

# forgatokonyvenkent egy minta a diagnosztikai abrakhoz
egy_minta <- function(seed = SEED) {
  set.seed(seed)
  mintak <- list()
  for (fk in FORGATOKONYVEK) {
    df <- fk$dgp(N)
    v  <- fk$psm[[1]]
    ps <- becsult_ps(df, v$valtozok)
    m  <- psm(df, ps, v$caliper)
    mintak[[fk$kod]] <- list(df = df, ps = ps, par = m$par, megtart = m$megtart)
  }
  mintak
}

main <- function() {
  dir.create(KIMENET, showWarnings = FALSE)

  mintak <- egy_minta()
  for (fk in FORGATOKONYVEK) {
    ment(abra_ps_eloszlas(fk, mintak[[fk$kod]]), fk, "1_ps_eloszlas")
    ment(abra_balance(fk, mintak[[fk$kod]]), fk, "2_balance", magassag = 3.6)
  }

  cat(sprintf("monte carlo: %d ismétlés, n = %d\n", ISMETLES, N))
  mc <- monte_carlo()
  write.csv(mc, file.path(KIMENET, "monte_carlo_nyers.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  for (fk in FORGATOKONYVEK) ment(abra_becslesek(fk, mc), fk, "3_becslesek")

  tabla <- osszefoglalo(mc)
  szamok <- sapply(tabla, is.numeric)
  tabla[szamok] <- round(tabla[szamok], 3)
  write.csv(tabla, file.path(KIMENET, "osszefoglalo.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  print(tabla[names(tabla) != "forgatokonyv"], row.names = FALSE)
  cat("kimenet:", KIMENET, "\n")
}

main()
