
`%||%` <- function(a, b) if (is.null(a)) b else a

with_seed <- function(seed, expr) {
  g <- globalenv(); had <- exists(".Random.seed", envir = g, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = g)
  on.exit(if (had) assign(".Random.seed", old, envir = g) else
    if (exists(".Random.seed", envir = g, inherits = FALSE)) rm(".Random.seed", envir = g))
  set.seed(seed); force(expr)
}

bh_by <- function(p, g, method = "BH") {
  out <- rep(NA_real_, length(p))
  for (u in unique(g)) { i <- which(g == u); out[i] <- stats::p.adjust(p[i], method) }
  out
}

GS_MIN_SIZE    <- 8
GS_MAX_SIZE    <- 500
GS_MIN_SIZE_SET <- c(
  KEGG_FATTY_ACID_BIOSYNTHESIS                                      = 5,
  REACTOME_NR1H2_NR1H3_REGULATE_GENE_EXPRESSION_LINKED_TO_LIPOGENESIS = 5)
TF_MIN_TARGETS <- 5
GS_MSIG_GMT_DIR <- NULL
TF_NETWORK_CSV  <- NULL

GS_KEYS <- c("tissue", "level", "subset", "comparison")
GS_MIN_DONORS_PER_GROUP <- if (exists("DEG_MIN_DONORS_PER_GROUP")) DEG_MIN_DONORS_PER_GROUP else 3
GS_MIN_DONORS_NUMERIC   <- 5

GS_MSIG_WANTED <- tibble::tribble(
  ~collection, ~subcollection, ~set_name,                                        ~panel,    ~counterpart,
  "H",  NA,            "HALLMARK_INTERFERON_GAMMA_RESPONSE",                     "immune",  NA,
  "H",  NA,            "HALLMARK_INTERFERON_ALPHA_RESPONSE",                     "immune",  NA,
  "H",  NA,            "HALLMARK_TNFA_SIGNALING_VIA_NFKB",                       "immune",  NA,
  "H",  NA,            "HALLMARK_INFLAMMATORY_RESPONSE",                         "immune",  NA,
  "H",  NA,            "HALLMARK_COMPLEMENT",                                    "immune",  NA,
  "H",  NA,            "HALLMARK_CHOLESTEROL_HOMEOSTASIS",                       "lipid",   NA,
  "H",  NA,            "HALLMARK_FATTY_ACID_METABOLISM",                         "lipid",   NA,
  "H",  NA,            "HALLMARK_PEROXISOME",                                    "lipid",   NA,
  "H",  NA,            "HALLMARK_ADIPOGENESIS",                                  "lipid",   NA,
  "H",  NA,            "HALLMARK_GLYCOLYSIS",                                    "metabolic", NA,
  "H",  NA,            "HALLMARK_MTORC1_SIGNALING",                              "metabolic", NA,
  "H",  NA,            "HALLMARK_HYPOXIA",                                       "metabolic", NA,
  "H",  NA,            "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",                     "stress",  NA,
  "H",  NA,            "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",               "stress",  NA,
  "H",  NA,            "HALLMARK_OXIDATIVE_PHOSPHORYLATION",                     "control", NA,
  "H",  NA,            "HALLMARK_MYC_TARGETS_V1",                                "control", NA,
  "H",  NA,            "HALLMARK_MYC_TARGETS_V2",                                "control", NA,
  "C2", "CP:REACTOME", "REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF", "sterol", NA,
  "C2", "CP:REACTOME", "REACTOME_PLASMA_LIPOPROTEIN_CLEARANCE",                  "lipid",   NA,
  "C2", "CP:REACTOME", "REACTOME_PEROXISOMAL_LIPID_METABOLISM",                  "lipid",   NA,
  "C2", "CP:REACTOME", "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",       "lipid",   NA,
  "C2", "CP:REACTOME", "REACTOME_SPHINGOLIPID_METABOLISM",                       "lipid",   NA,
  "C2", "CP:REACTOME", "REACTOME_BETA_OXIDATION_OF_VERY_LONG_CHAIN_FATTY_ACIDS",  "VLCFA",   "PEROXISOMAL_BETA_OXIDATION",
  "C2", "CP:KEGG_LEGACY", "KEGG_LYSOSOME",                                       "lipid",   "LYSOSOMAL_DEGRADATIVE",
  "C2", "CP:REACTOME", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",             "control", NA,
  "C2", "CP:REACTOME",    "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP", "sterol", "SREBP2_TARGETS_CORE",
  "C2", "CP:KEGG_LEGACY", "KEGG_TERPENOID_BACKBONE_BIOSYNTHESIS",                  "sterol", "MEVALONATE_PRE_SQUALENE",
  "C2", "CP:REACTOME",    "REACTOME_CHOLESTEROL_BIOSYNTHESIS",                     "sterol", "STEROL_POST_SQUALENE",
  "C2", "CP:REACTOME",    "REACTOME_NR1H2_NR1H3_REGULATE_GENE_EXPRESSION_LINKED_TO_LIPOGENESIS", "lipid", "SREBF1_LIPOGENIC",
  "C2", "CP:KEGG_LEGACY", "KEGG_FATTY_ACID_BIOSYNTHESIS",                          "lipid",  "SREBF1_LIPOGENIC",
  "C2", "CP:REACTOME",    "REACTOME_NR1H3_NR1H2_REGULATE_GENE_EXPRESSION_LINKED_TO_CHOLESTEROL_TRANSPORT_AND_EFFLUX", "lipid", "LXR_TARGETS_STRICT",
  "C2", "CP:REACTOME",    "REACTOME_THE_NLRP3_INFLAMMASOME",                       "stress", "INFLAMMASOME_PRIMING"
)

GS_REGULATOR_GENES <- c(
  "NR1H2", "NR1H3", "RXRA", "RXRB", "RXRG", "NRIP1", "PPARA",
  "SREBF1", "SREBF2", "NFYA", "NFYB", "NFYC", "SP1", "MTF1",
  "NCOA1", "NCOA2", "NCOA3", "NCOA6", "NCOR1", "NCOR2", "EP300", "CREBBP",
  "CARM1", "MED1", "SMARCD3", "CHD9", "HELZ2", "TGS1", "TBL1X", "TBL1XR1",
  "HDAC3", "GPS2", "KDM1A", "KDM1B", "KDM3A", "KDM4A")
GS_MSIG_EXCLUDE <- list(
  REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP = GS_REGULATOR_GENES,
  REACTOME_NR1H3_NR1H2_REGULATE_GENE_EXPRESSION_LINKED_TO_CHOLESTEROL_TRANSPORT_AND_EFFLUX = GS_REGULATOR_GENES)

GS_CUSTOM <- list(
  SREBP2_TARGETS_CORE = c(
    "HMGCS1", "HMGCR", "MVK", "PMVK", "MVD", "IDI1", "FDPS", "FDFT1",
    "SQLE", "LSS", "CYP51A1", "DHCR7", "DHCR24", "ACAT2",
    "LDLR", "PCSK9", "INSIG1", "STARD4", "NPC1", "ACSS2"),
  MEVALONATE_PRE_SQUALENE = c(
    "ACAT2", "HMGCS1", "HMGCR", "MVK", "PMVK", "MVD",
    "IDI1", "FDPS", "GGPS1", "FDFT1"),
  STEROL_POST_SQUALENE = c(
    "SQLE", "LSS", "CYP51A1", "MSMO1", "NSDHL", "HSD17B7", "SC5D",
    "EBP", "DHCR7", "DHCR24"),

  SREBF1_LIPOGENIC = c(
    "FASN", "ACACA", "ACLY", "SCD", "SCD5", "ELOVL6", "ELOVL5",
    "ACSS2", "AGPAT2", "ACSL3", "DGAT1", "DGAT2", "SREBF1"),

  PEROXISOMAL_BETA_OXIDATION = c(
    "ABCD3", "ACOX1", "ACOX3", "ACAA1", "HSD17B4", "SCP2", "CROT", "ACOT8"),

  LAM_CORE = c(
    "TREM2", "TYROBP", "GPNMB", "SPP1", "CD9", "LGALS3", "CD36",
    "LPL", "FABP5", "APOE", "APOC1", "ITGAX", "CST7",
    "CLEC7A", "AXL", "MSR1"),
  LYSOSOMAL_DEGRADATIVE = c(
    "CTSB", "CTSD", "CTSL", "CTSZ", "CTSS", "PSAP", "GRN", "HEXA", "HEXB",
    "ASAH1", "GLA", "GBA", "CD63", "ACP5", "LAMP1", "LAMP2", "TFEB"),

  LXR_TARGETS = c(
    "ABCA1", "ABCG1", "APOE", "APOC1", "MYLIP", "PLTP", "ARL4C",
    "NR1H3", "SREBF1", "SCD", "LPL"),
  LXR_TARGETS_STRICT = c(
    "ABCA1", "ABCG1", "MYLIP", "PLTP", "ARL4C"),

  INFLAMMASOME_PRIMING = c(
    "NLRP3", "PYCARD", "CASP1", "GSDMD", "IL1B", "IL18", "NEK7",
    "CARD8", "P2RX7"),

  ER_STRESS_UPR_TARGETS = c(
    "HSPA5", "DDIT3", "HERPUD1", "DNAJB9", "SDF2L1", "SEL1L", "HYOU1",
    "EDEM1", "DNAJC3", "TRIB3", "WFS1", "CHAC1"),

  IEG_DISSOCIATION_STRESS = c(
    "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
    "ATF3", "IER2", "IER3", "KLF4", "KLF6", "ZFP36",
    "BTG1", "BTG2", "PPP1R15A", "RHOB", "SOCS3", "NR4A1",
    "HSPA1A", "HSPA1B", "MCL1"),

  GC_RESPONSE_CANONICAL = c(
    "FKBP5", "TSC22D3", "ZBTB16", "KLF9", "PER1", "DDIT4", "SESN1",
    "TXNIP", "IRS2", "THBD", "SLA", "CEBPD", "ADORA3", "PIK3R1",
    "SAP30", "MT1X", "AMPD3"),
  GC_RESPONSE_MYELOID = c(
    "IL1R2", "FPR1", "ADORA3", "AMPD3", "ALOX15B", "SLC1A3",
    "FKBP5", "TSC22D3", "ZBTB16"),

  EFFEROCYTOSIS_RECEPTORS = c(
    "MERTK", "AXL", "TYRO3", "ELMO1", "DOCK1", "ITGB5", "ITGAV",
    "CD36", "MARCO"),
  MHCII_ANTIGEN_PRESENTATION = c(
    "HLA-DRA", "HLA-DMA", "HLA-DMB", "CD74", "CIITA",
    "CTSS", "LGMN", "IFI30", "RFX5", "RFXAP")
)
GS_CUSTOM_PANEL <- c(
  SREBP2_TARGETS_CORE = "sterol", MEVALONATE_PRE_SQUALENE = "sterol",
  STEROL_POST_SQUALENE = "sterol",
  SREBF1_LIPOGENIC = "lipid", LAM_CORE = "lipid", LYSOSOMAL_DEGRADATIVE = "lipid",
  LXR_TARGETS = "lipid", LXR_TARGETS_STRICT = "lipid",
  PEROXISOMAL_BETA_OXIDATION = "VLCFA",
  INFLAMMASOME_PRIMING = "stress", ER_STRESS_UPR_TARGETS = "stress",
  EFFEROCYTOSIS_RECEPTORS = "immune", MHCII_ANTIGEN_PRESENTATION = "immune",
  GC_RESPONSE_CANONICAL = "confound", GC_RESPONSE_MYELOID = "confound",
  IEG_DISSOCIATION_STRESS = "control")

GS_CUSTOM_RETIRED <- names(GS_CUSTOM)
GS_GC_GENES <- c("FKBP5", "TSC22D3", "ZBTB16", "DDIT4", "PER1")

GS_PUBLISHED_DIR <- here::here("genesets")
gs_published_sets <- function(dir = GS_PUBLISHED_DIR) {
  fs <- list.files(dir, "\\.txt$", full.names = TRUE)
  out <- list(); panel <- character(0); src <- character(0)
  for (f in fs) {
    l   <- trimws(readLines(f, warn = FALSE))
    hdr <- grep("^#", l, value = TRUE)
    val <- function(k) sub(paste0("^#\\s*", k, ":\\s*"), "", grep(paste0("^#\\s*", k, ":"), hdr, value = TRUE))[1]
    g   <- unique(l[nzchar(l) & !grepl("^#", l)])
    nm  <- toupper(tools::file_path_sans_ext(basename(f)))
    if (!length(g)) { message("published gene set ", nm, ": no genes in ", basename(f), " -- not run"); next }
    if (identical(tolower(val("species")), "mouse")) {
      if (!requireNamespace("babelgene", quietly = TRUE)) stop("babelgene needed for mouse gene set ", nm, call. = FALSE)
      o <- babelgene::orthologs(g, species = "mouse", human = FALSE)
      miss <- setdiff(g, o$symbol)
      if (length(miss)) message(nm, ": no human ortholog for ", toString(miss))
      g <- unique(o$human_symbol)
    }
    out[[nm]] <- g
    panel[[nm]] <- if (is.na(val("panel"))) "published" else val("panel")
    src[[nm]] <- val("source")
  }
  attr(out, "panel") <- panel; attr(out, "source") <- src
  out
}
gs_custom_active <- function() {
  pub <- gs_published_sets()
  c(GS_CUSTOM[setdiff(names(GS_CUSTOM), GS_CUSTOM_RETIRED)], pub)
}
GS_CUSTOM_PANEL <- c(GS_CUSTOM_PANEL, KEREN_SHAUL_DAM = "lipid", VAN_DEN_BRINK_DISSOCIATION = "control")
local({
  miss <- setdiff(names(GS_CUSTOM), names(GS_CUSTOM_PANEL))
  if (length(miss)) warning("GS_CUSTOM sets without a panel: ", toString(miss), call. = FALSE)
  bad <- setdiff(stats::na.omit(GS_MSIG_WANTED$counterpart), names(GS_CUSTOM))
  if (length(bad)) warning("GS_MSIG_WANTED counterparts not in GS_CUSTOM: ", toString(bad),
                           call. = FALSE)
})

gs_control_names <- function()
  unique(c("RIBOSOME_RP", "OXPHOS_SUBUNITS", "MITO_RIBOSOME",
           names(GS_CUSTOM_PANEL)[GS_CUSTOM_PANEL == "control"],
           GS_MSIG_WANTED$set_name[GS_MSIG_WANTED$panel == "control"]))

GS_STRIP_TECHNICAL <- FALSE
GS_RUN_STRIP_SENS  <- TRUE

GS_CENTRE_STAT     <- FALSE
GS_CENTRE_SENS     <- TRUE
GS_RESID_SENS      <- TRUE
GS_GSEA_PARAM_SENS <- 0.5
GS_FGSEA_EPS       <- 1e-10
GS_CAMERA_COR      <- 0.01
GS_CAMERA_COR_SENS <- 0.05
GS_MATCH_N    <- 1000
GS_MATCH_BINS <- 20

GS_VOOM_BLOCK <- TRUE

GS_IUT_COLS <- c("pval_camera_vc_est", "pval_fgsea")
GS_DIR_OF   <- c(pval_fgsea = "NES", pval_camera = "camera_direction",
                 pval_camera_vc = "direction_camera_vc", pval_camera_vc_est = "direction_camera_vc_est",
                 pval_fry = "direction_fry")

TF_OF_INTEREST <- c("SREBF1","SREBF2","NR1H2","NR1H3","PPARG","PPARA","NR3C1",
                    "STAT1","IRF1","NFKB1","RELA","CEBPB","TFEB","MITF","MAFB",
                    "HIF1A","ATF4","XBP1","NFE2L2")

gs_read_gmt <- function(f) {
  l <- strsplit(readLines(f, warn = FALSE), "\t")
  stats::setNames(lapply(l, function(x) unique(x[-(1:2)])), vapply(l, `[`, "", 1))
}

gs_fetch_msigdb <- function(wanted = GS_MSIG_WANTED, gmt_dir = GS_MSIG_GMT_DIR,
                            cache = file.path(path$objects, "msigdb_selected_sets.rds")) {
  if (!is.null(gmt_dir) && dir.exists(gmt_dir)) {
    raw  <- unlist(lapply(list.files(gmt_dir, "\\.gmt$", full.names = TRUE), gs_read_gmt),
                   recursive = FALSE)
    sets <- raw[intersect(names(raw), wanted$set_name)]
    src  <- paste("local .gmt:", gmt_dir)
  } else {
    sets <- if (file.exists(cache)) tryCatch(readRDS(cache), error = function(e) list()) else list()
    src  <- paste("cache:", cache)
    if (!all(wanted$set_name %in% names(sets)) && requireNamespace("msigdbr", quietly = TRUE)) {
      get1 <- function(coll, sub) {
        subs <- if (is.na(sub)) list(NULL) else
          as.list(unique(c(sub, if (sub == "CP:KEGG_LEGACY") "CP:KEGG")))
        for (s in subs) {
          a <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = coll,
                                         subcollection = s), error = function(e) NULL)
          if (is.null(a) || !nrow(a))
            a <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = coll,
                                           subcategory = s), error = function(e) NULL)
          if (!is.null(a) && nrow(a)) return(a[, c("gs_name", "gene_symbol")])
        }
        message("msigdbr: nothing returned for ", coll, " / ", ifelse(is.na(sub), "-", sub))
        NULL
      }
      combos <- unique(wanted[, c("collection", "subcollection")])
      df <- dplyr::bind_rows(lapply(seq_len(nrow(combos)), function(i)
        get1(combos$collection[i], combos$subcollection[i])))
      if (nrow(df)) {
        df   <- df[df$gs_name %in% wanted$set_name, ]
        new  <- lapply(split(df$gene_symbol, df$gs_name), unique)
        sets <- c(sets[setdiff(names(sets), names(new))], new)
        src  <- paste("msigdbr", utils::packageVersion("msigdbr"))
        tryCatch(saveRDS(sets, cache), error = function(e)
          warning("could not cache MSigDB sets: ", conditionMessage(e), call. = FALSE))
      }
    }
  }
  sets <- sets[intersect(wanted$set_name, names(sets))]
  miss <- setdiff(wanted$set_name, names(sets))
  if (length(miss)) message("MSigDB sets not available (", src, "): ", toString(miss))
  attr(sets, "source") <- src
  sets
}

gs_control_sets <- function(universe) list(
  RIBOSOME_RP     = if (exists("DEG_RIBO_GENES")) intersect(DEG_RIBO_GENES, universe) else
    grep("^RP[LS][0-9]+[AXY]?$|^RPLP[0-9]$|^RPSA$", universe, value = TRUE),
  OXPHOS_SUBUNITS = grep("^NDUF|^COX[0-9]|^UQCR|^ATP5|^SDH[ABCD]$|^CYC1$", universe, value = TRUE),
  MITO_RIBOSOME   = grep(if (exists("DEG_MITO_RIBO_PATTERN")) DEG_MITO_RIBO_PATTERN else "^MRP[LS]",
                         universe, value = TRUE))

strip_technical <- function(sets, tech_sets) {
  tech <- unique(unlist(tech_sets, use.names = FALSE))
  lapply(sets, setdiff, tech)
}

gs_exclude_regulators <- function(sets, exclude = GS_MSIG_EXCLUDE) {
  n <- stats::setNames(rep(0L, length(sets)), names(sets))
  for (s in intersect(names(exclude), names(sets))) {
    before <- length(sets[[s]])
    sets[[s]] <- setdiff(sets[[s]], exclude[[s]])
    n[s] <- before - length(sets[[s]])
  }
  attr(sets, "n_removed") <- n
  sets
}

gs_assemble <- function(universe, msig = gs_fetch_msigdb(), strip = GS_STRIP_TECHNICAL) {
  src   <- attr(msig, "source") %||% "MSigDB"
  msig  <- gs_exclude_regulators(msig)
  n_reg <- attr(msig, "n_removed")
  ctrl  <- gs_control_sets(universe)
  cust  <- gs_custom_active()
  other <- c(msig, cust)
  n_before <- lengths(other)
  if (strip) {
    panel_of <- c(GS_MSIG_WANTED$panel[match(names(msig), GS_MSIG_WANTED$set_name)],
                  unname(GS_CUSTOM_PANEL[names(cust)]))
    bio <- !panel_of %in% "control"
    other[bio] <- strip_technical(other[bio], ctrl)
  }
  sets <- c(other, ctrl)
  mi   <- match(names(msig), GS_MSIG_WANTED$set_name)
  prov <- dplyr::bind_rows(
    tibble::tibble(set = names(msig), source = src,
                   panel = GS_MSIG_WANTED$panel[mi],
                   counterpart = GS_MSIG_WANTED$counterpart[mi]),
    tibble::tibble(set = as.character(names(cust)),
                   source = as.character(ifelse(names(cust) %in% names(GS_CUSTOM), "curated",
                                   paste("published:", attr(gs_published_sets(), "source")[names(cust)]))),
                   panel = as.character(unname(GS_CUSTOM_PANEL[names(cust)]))),
    tibble::tibble(set = names(ctrl), source = "regex over tested genes", panel = "control"))
  prov$n_genes_set    <- lengths(sets)[prov$set]
  prov$n_genes_tested <- vapply(sets[prov$set], function(g) length(intersect(g, universe)), 0L)
  prov$n_regulators_removed <- unname(ifelse(prov$set %in% names(n_reg), n_reg[prov$set],
                                             NA_integer_))
  prov$n_technical_removed <- unname(ifelse(prov$set %in% names(n_before),
                                            n_before[prov$set] - lengths(sets)[prov$set], NA_integer_))
  list(sets = sets, provenance = prov)
}

gs_rank_vector <- function(df, centre = GS_CENTRE_STAT, residualise = FALSE) {
  df <- df[!is.na(df$stat) & !is.na(df$gene), , drop = FALSE]
  nd <- sum(duplicated(df$gene))
  if (nd) warning(nd, " duplicated gene symbols; largest |stat| kept -- fix the ID mapping",
                  call. = FALSE)
  df <- df[order(-abs(df$stat)), , drop = FALSE]
  df <- df[!duplicated(df$gene), , drop = FALSE]
  shift <- stats::median(df$stat)
  val <- if (centre) df$stat - shift else df$stat
  if (residualise) {
    if (!"baseMean" %in% names(df)) stop("residualise = TRUE needs baseMean", call. = FALSE)
    x   <- log10(df$baseMean + 1)
    fit <- stats::loess(df$stat ~ x, span = 0.3, family = "symmetric", degree = 1)
    val <- df$stat - stats::predict(fit, x)
  }
  out <- stats::setNames(val, df$gene)
  attr(out, "shift") <- shift
  if ("baseMean" %in% names(df)) attr(out, "baseMean") <- stats::setNames(df$baseMean, df$gene)
  out
}

gs_matched_null <- function(stats, sets, n = GS_MATCH_N, bins = GS_MATCH_BINS, seed = 1) {
  bm <- attr(stats, "baseMean")
  if (is.null(bm)) return(stats::setNames(rep(NA_real_, length(sets)), names(sets)))
  bm  <- bm[names(stats)]
  bin <- stats::setNames(cut(rank(bm, ties.method = "first"), bins, labels = FALSE), names(stats))
  pools <- split(names(stats), bin)
  with_seed(seed, vapply(sets, function(g) {
    g <- intersect(g, names(stats))
    if (!length(g)) return(NA_real_)
    obs <- mean(stats[g])
    k   <- table(factor(bin[g], levels = names(pools)))
    tot <- numeric(n)
    for (b in names(k)[k > 0]) {
      p <- stats[setdiff(pools[[b]], g)]
      if (!length(p)) p <- stats[pools[[b]]]
      tot <- tot + rowSums(matrix(p[sample.int(length(p), n * k[[b]], replace = TRUE)], nrow = n))
    }
    null <- tot / length(g)
    m0   <- mean(null)
    (1 + sum(abs(null - m0) >= abs(obs - m0))) / (1 + n)
  }, numeric(1)))
}

gs_test_family <- function(stats, sets, seed = 1, control_sets = gs_control_names(),
                           gsea_param_sens = GS_GSEA_PARAM_SENS, matched = TRUE) {
  sets  <- lapply(sets, intersect, names(stats))
  min_k <- ifelse(names(sets) %in% names(GS_MIN_SIZE_SET),
                  GS_MIN_SIZE_SET[names(sets)], GS_MIN_SIZE)
  sets  <- sets[lengths(sets) >= min_k & lengths(sets) <= GS_MAX_SIZE]
  if (!length(sets)) return(NULL)
  out <- tibble::tibble(set = names(sets), size = lengths(sets),
                        is_control = names(sets) %in% control_sets,
                        NES = NA_real_, pval_fgsea = NA_real_, padj_fgsea = NA_real_,
                        NES_sens = NA_real_, pval_fgsea_sens = NA_real_, padj_fgsea_sens = NA_real_,
                        leading_edge = NA_character_)
  st <- stats::setNames(as.numeric(stats), names(stats))
  if (requireNamespace("fgsea", quietly = TRUE)) {
    fg1 <- function(gp) as.data.frame(suppressWarnings(fgsea::fgseaMultilevel(
      sets, st, minSize = min(lengths(sets)), maxSize = GS_MAX_SIZE,
      eps = GS_FGSEA_EPS, gseaParam = gp, BPPARAM = BiocParallel::SerialParam())))
    fg <- with_seed(seed, fg1(1))
    m  <- match(out$set, fg$pathway)
    out$NES <- fg$NES[m]; out$pval_fgsea <- fg$pval[m]
    out$leading_edge <- vapply(fg$leadingEdge[m], function(x)
      if (is.null(x)) NA_character_ else paste(utils::head(x, 15), collapse = ","), character(1))
    if (!is.null(gsea_param_sens)) {
      f2 <- with_seed(seed, fg1(gsea_param_sens)); m2 <- match(out$set, f2$pathway)
      out$NES_sens <- f2$NES[m2]; out$pval_fgsea_sens <- f2$pval[m2]
    }
  }
  idx <- lapply(sets, function(s) match(s, names(st)))
  cm  <- limma::cameraPR(st, idx)
  cs  <- limma::cameraPR(st, idx, inter.gene.cor = GS_CAMERA_COR_SENS)
  out$camera_direction <- cm[out$set, "Direction"]
  out$pval_camera      <- cm[out$set, "PValue"]
  out$pval_camera_sens <- cs[out$set, "PValue"]
  out$pval_matched     <- if (matched) unname(gs_matched_null(stats, sets[out$set], seed = seed)) else NA_real_
  for (p in c("fgsea", "fgsea_sens", "camera", "camera_sens", "matched"))
    out[[paste0("padj_", p)]] <- bh_by(out[[paste0("pval_", p)]], out$is_control)
  out$stat_shift <- attr(stats, "shift") %||% NA_real_
  out
}

.gs_voom_prep <- function(counts, cd, design, group = NULL) {
  des <- stats::model.matrix(design, cd)
  if (nrow(des) - ncol(des) < 2) stop("too few residual df for voom", call. = FALSE)
  dge <- edgeR::DGEList(counts)
  dge <- dge[edgeR::filterByExpr(dge, design = des), , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  g   <- if (isTRUE(GS_VOOM_BLOCK) && !is.null(group)) droplevels(factor(cd[[group]])) else NULL
  if (!is.null(g) && min(table(g)) < 2) g <- NULL
  list(v = limma::voomWithQualityWeights(dge, des, var.group = g, plot = FALSE), des = des)
}

.gs_contrast_vector <- function(cd, design, ct) {
  if (length(ct) == 1) {
    mm <- stats::model.matrix(design, cd)
    if (!ct %in% colnames(mm)) stop("coefficient '", ct, "' not in the design", call. = FALSE)
    return(stats::setNames(as.numeric(colnames(mm) == ct), colnames(mm)))
  }
  tt   <- stats::delete.response(stats::terms(design))
  vars <- intersect(all.vars(tt), names(cd))
  xl   <- lapply(Filter(is.factor, cd[vars]), levels)
  nd   <- cd[c(1, 1), , drop = FALSE]
  nd[[ct[1]]] <- factor(c(ct[2], ct[3]), levels = levels(cd[[ct[1]]]))
  mm <- stats::model.matrix(tt, stats::model.frame(tt, nd, xlev = xl))
  mm[1, ] - mm[2, ]
}

.gs_voom_core <- function(counts, cd, design, ct, sets) {
  p   <- .gs_voom_prep(counts, cd, design, group = if (length(ct) == 3) ct[1] else NULL)
  cv  <- .gs_contrast_vector(cd, design, ct)
  idx <- limma::ids2indices(sets, rownames(p$v), remove.empty = FALSE)
  min_k <- ifelse(names(idx) %in% names(GS_MIN_SIZE_SET), GS_MIN_SIZE_SET[names(idx)], GS_MIN_SIZE)
  idx <- idx[lengths(idx) >= min_k & lengths(idx) <= GS_MAX_SIZE]
  if (!length(idx)) return(NULL)
  cam <- limma::camera(p$v, idx, p$des, contrast = cv, inter.gene.cor = GS_CAMERA_COR)
  ce  <- limma::camera(p$v, idx, p$des, contrast = cv, inter.gene.cor = NA)
  fr  <- limma::fry(p$v, idx, p$des, contrast = cv)
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(p$v, p$des), matrix(cv, ncol = 1)))
  b <- fit$coefficients[, 1]; tt <- fit$t[, 1]
  nm <- names(idx)
  tibble::tibble(set = nm, size = lengths(idx),
                 mean_logFC = vapply(idx, function(i) mean(b[i]), 0),
                 mean_t     = vapply(idx, function(i) mean(tt[i]), 0),
                 direction_camera_vc = cam[nm, "Direction"], pval_camera_vc = cam[nm, "PValue"],
                 corr_camera_vc_est = ce[nm, "Correlation"],
                 direction_camera_vc_est = ce[nm, "Direction"], pval_camera_vc_est = ce[nm, "PValue"],
                 direction_fry = fr[nm, "Direction"], pval_fry = fr[nm, "PValue"])
}

.gs_check_groups <- function(cd, ct) {
  if (length(ct) == 1) {
    if (nrow(cd) < GS_MIN_DONORS_NUMERIC)
      stop(sprintf("< %d donors for the numeric coefficient %s", GS_MIN_DONORS_NUMERIC, ct), call. = FALSE)
    return(invisible(TRUE))
  }
  tab <- table(cd[[ct[1]]])[ct[2:3]]
  if (any(is.na(tab) | tab < GS_MIN_DONORS_PER_GROUP))
    stop(sprintf("group(s) with < %d donors: %s", GS_MIN_DONORS_PER_GROUP,
                 paste(names(tab), tab, sep = "=", collapse = ", ")), call. = FALSE)
}

gs_test_family_voom <- function(dds, ct, sets, control_sets = gs_control_names())
  .gs_test_voom(DESeq2::counts(dds), as.data.frame(SummarizedExperiment::colData(dds)),
                DESeq2::design(dds), ct, sets, control_sets)

.gs_test_voom <- function(counts, cd, design, ct, sets, control_sets = gs_control_names()) {
  .gs_check_groups(cd, ct)
  r <- .gs_voom_core(counts, cd, design, ct, sets)
  if (is.null(r)) return(NULL)
  r$is_control     <- r$set %in% control_sets
  r$padj_camera_vc <- bh_by(r$pval_camera_vc, r$is_control)
  r$padj_camera_vc_est <- bh_by(r$pval_camera_vc_est, r$is_control)
  r$padj_fry       <- bh_by(r$pval_fry, r$is_control)
  r
}

gs_registry_add <- function(reg, tissue, level, subset, dds, contrasts, qc = NULL)
  c(reg, list(list(tissue = tissue, level = level, subset = subset, dds = dds,
                   contrasts = contrasts, qc = qc)))

gs_run_voom <- function(registry, sets, comparisons = NULL) {
  dplyr::bind_rows(lapply(registry, function(e)
    dplyr::bind_rows(lapply(names(e$contrasts), function(cn) {
      if (!is.null(comparisons) && !cn %in% comparisons) return(NULL)
      r <- tryCatch(gs_test_family_voom(e$dds, e$contrasts[[cn]], sets),
                    error = function(err) {
                      message("gs_run_voom skipped ", e$tissue, " | ", e$level, " | ", e$subset,
                              " | ", cn, ": ", conditionMessage(err)); NULL })
      if (is.null(r)) return(NULL)
      k <- tibble::tibble(tissue = e$tissue, level = e$level, subset = e$subset, comparison = cn)
      dplyr::bind_cols(k[rep(1, nrow(r)), ], r)
    }))))
}

gs_run <- function(deg_tbl, sets, comparisons, levels_use = c("all_myeloid", "celltype"),
                   centre_sens = GS_CENTRE_SENS, resid_sens = GS_RESID_SENS, light = FALSE) {
  if (!requireNamespace("fgsea", quietly = TRUE))
    warning("fgsea not installed -- cameraPR only. BiocManager::install('fgsea')", call. = FALSE)
  fam_rep <- deg_family_report(deg_tbl, comparisons, levels_use, "gs_run")
  fam <- fam_rep %>% dplyr::filter(is.na(reason)) %>%
    dplyr::select(tissue, level, subset, comparison)
  variants <- c("raw", if (centre_sens && !light) "median_centred",
                if (resid_sens && !light) "intensity_resid")
  dplyr::bind_rows(lapply(seq_len(nrow(fam)), function(i) {
    f  <- fam[i, ]
    df <- deg_tbl %>% dplyr::filter(status == "ok", tissue == f$tissue, level == f$level,
                                    subset == f$subset, comparison == f$comparison)
    dplyr::bind_rows(lapply(variants, function(v) {
      rv <- tryCatch(gs_rank_vector(df, centre = v == "median_centred",
                                    residualise = v == "intensity_resid"),
                     error = function(e) { message("gs_run ", v, ": ", conditionMessage(e)); NULL })
      if (is.null(rv)) return(NULL)
      r <- gs_test_family(rv, sets, matched = v == "raw" && !light,
                          gsea_param_sens = if (v == "raw" && !light) GS_GSEA_PARAM_SENS else NULL)
      if (is.null(r)) return(NULL)
      r$ranking <- v
      dplyr::bind_cols(f[rep(1, nrow(r)), ], r)
    }))
  }))
}

.gs_dir <- function(x) if (is.numeric(x)) sign(x) else
  ifelse(x == "Up", 1, ifelse(x == "Down", -1, NA_real_))

gs_iut <- function(res, p_cols = GS_IUT_COLS) {
  need <- c(p_cols, unname(GS_DIR_OF[p_cols]), "is_control")
  if (length(miss <- setdiff(need, names(res)))) stop("gs_iut: missing ", toString(miss), call. = FALSE)
  P <- as.matrix(res[p_cols])
  D <- do.call(cbind, lapply(p_cols, function(p) .gs_dir(res[[GS_DIR_OF[[p]]]])))
  agree <- apply(D, 1, function(d) length(unique(d[!is.na(d)])) <= 1)
  p <- apply(P, 1, max)
  p[!agree & !is.na(p)] <- 1
  res$p_iut <- p
  keys <- intersect(c(GS_KEYS, "ranking"), names(res))
  grp  <- do.call(paste, c(lapply(c(keys, "is_control"), function(k) res[[k]]), sep = "|"))
  res$padj_iut <- bh_by(res$p_iut, grp)
  res
}

gs_combine <- function(res_gene, res_voom, iut_cols = GS_IUT_COLS) {
  g <- res_gene
  if ("ranking" %in% names(g)) { g <- g[g$ranking == "raw", ]; g$ranking <- NULL }
  g <- g[, setdiff(names(g), "is_control")]
  names(g)[names(g) == "size"] <- "size_gene"
  gs_iut(dplyr::inner_join(res_voom, g, by = c(GS_KEYS, "set")), iut_cols)
}

gs_global_padj <- function(res, p_col = "p_iut", by = "comparison", method = "BH") {
  res$padj_global <- bh_by(res[[p_col]], paste(res[[by]], res$is_control, sep = "|"), method)
  res
}

gs_independent <- function(res, sets, p_col = "p_iut", alpha = 0.05, overlap = 0.5) {
  jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  res$independent <- NA
  fam <- do.call(paste, c(lapply(intersect(c(GS_KEYS, "ranking"), names(res)),
                                 function(k) res[[k]]), sep = "|"))
  for (f in unique(fam)) {
    i <- which(fam == f & !is.na(res[[p_col]]) & res[[p_col]] < alpha & !res$is_control)
    i <- i[order(res[[p_col]][i])]
    kept <- character(0)
    for (j in i) {
      s <- res$set[j]
      ok <- !length(kept) || all(vapply(kept, function(k) jac(sets[[s]], sets[[k]]), 0) < overlap)
      res$independent[j] <- ok
      if (ok) kept <- c(kept, s)
    }
  }
  res
}

gs_agree <- function(res, padj_max = 0.05) {
  if ("agree" %in% names(res)) return(res$agree %in% TRUE)
  warning("no agree column: legacy fgsea+cameraPR rule; use gs_finalise()", call. = FALSE)
  !is.na(res$padj_fgsea) & res$padj_fgsea < padj_max & res$padj_camera < padj_max &
    sign(res$NES) == ifelse(res$camera_direction == "Up", 1, -1)
}

gs_shift_robust <- function(res, padj_max = 0.05, na_pass = FALSE) {
  chk <- function(p) if (na_pass) is.na(p) | p < padj_max else !is.na(p) & p < padj_max
  gs_agree(res, padj_max) & chk(res$padj_camera_sens) & chk(res$padj_matched)
}

gs_counterpart_table <- function(res, provenance, sets = NULL) {
  pairs <- provenance[!is.na(provenance$counterpart), c("set", "counterpart")]
  if (!nrow(pairs)) return(tibble::tibble())
  keys <- GS_KEYS
  pub  <- dplyr::inner_join(res, pairs, by = "set")
  cur  <- res[res$set %in% pairs$counterpart, c(keys, "set", "NES", "padj_camera")]
  names(cur)[names(cur) %in% c("set", "NES", "padj_camera")] <-
    c("counterpart", "NES_curated", "padj_camera_curated")
  out <- dplyr::left_join(pub[, c(keys, "set", "counterpart", "NES", "padj_camera")], cur,
                          by = c(keys, "counterpart")) %>%
    dplyr::rename(published = set, NES_published = NES,
                  padj_camera_published = padj_camera) %>%
    dplyr::mutate(same_direction = sign(NES_published) == sign(NES_curated))
  if (!is.null(sets))
    out$jaccard_genes <- mapply(function(a, b) {
      x <- sets[[a]]; y <- sets[[b]]
      if (is.null(x) || is.null(y)) NA_real_ else length(intersect(x, y)) / length(union(x, y))
    }, out$published, out$counterpart)
  out
}

gs_permutation_null <- function(pb, sets, mode = c("controls_split", "permute"), n_perm = 200,
                                seed = 1, sf_mode = DEG_SF_MODE, alpha = 0.05) {
  mode <- match.arg(mode)
  if (!exists("deg_null_labelings") || !"sf_mode" %in% names(formals(fit_dds)))
    stop("gs_permutation_null() needs R/deg_fixes.R (deg_null_labelings(), fit_dds(sf_mode)) ",
         "sourced after R/03_deg.R", call. = FALSE)
  sm <- pb$meta
  sm$diagnosis <- factor(sm$diagnosis, levels = c("control", "X-ALD"))
  sm <- sm[!is.na(sm$diagnosis), , drop = FALSE]
  n_ctrl <- sum(sm$diagnosis == "control")
  if (mode == "controls_split" && n_ctrl < 2 * GS_MIN_DONORS_PER_GROUP)
    stop(sprintf("controls_split: %d controls give %d vs %d, below the floor of %d per group -- use mode = 'permute'",
                 n_ctrl, floor(n_ctrl / 2), ceiling(n_ctrl / 2), GS_MIN_DONORS_PER_GROUP), call. = FALSE)
  labs <- deg_null_labelings(sm, mode, n_perm, seed)
  true_lab <- stats::setNames(as.character(sm$diagnosis), rownames(sm))
  n_swapped <- vapply(labs, function(m) sum(as.character(m$diagnosis) != true_lab[rownames(m)], na.rm = TRUE), 0)
  ct   <- c("diagnosis", "X-ALD", "control")
  pc   <- c("pval_fgsea", "pval_camera", "pval_matched", "pval_camera_vc", "pval_camera_vc_est", "pval_fry")
  one <- lapply(seq_along(labs), function(j) {
    m <- labs[[j]]
    tryCatch({
      check_min_donors(m, "diagnosis")
      dds <- fit_dds(pb$counts[, rownames(m), drop = FALSE], m, ~ diagnosis, sf_mode)
      de  <- get_res(dds, ct, shrink = "none")
      g <- gs_test_family(gs_rank_vector(data.frame(gene = rownames(de), stat = de$stat,
                                                    baseMean = de$baseMean), centre = FALSE), sets)
      v <- gs_test_family_voom(dds, ct, sets)
      out <- dplyr::full_join(g[, c("set", "is_control", intersect(pc, names(g)))],
                              v[, c("set", intersect(pc, names(v)))], by = "set")
      out$perm <- j; out$n_swapped <- n_swapped[j]; out
    }, error = function(e) { message("perm ", j, " failed: ", conditionMessage(e)); NULL })
  })
  all  <- dplyr::bind_rows(one)
  nc   <- all[!all$is_control %in% TRUE, ]
  bal_min <- ceiling(min(table(sm$diagnosis)) / 2)
  nc_bal  <- nc[nc$n_swapped >= bal_min, ]
  cols <- intersect(pc, names(nc))
  summ <- dplyr::bind_rows(lapply(cols, function(p) {
    any_sig <- tapply(nc[[p]], nc$perm, function(x) any(stats::p.adjust(x, "BH") < alpha, na.rm = TRUE))
    tibble::tibble(test = sub("^pval_", "", p),
                   rate_p_lt_alpha = mean(nc[[p]] < alpha, na.rm = TRUE),
                   rate_p_lt_alpha_balanced = mean(nc_bal[[p]] < alpha, na.rm = TRUE),
                   frac_perms_any_padj = mean(any_sig),
                   n_perm = length(any_sig), n_perm_balanced = length(unique(nc_bal$perm)))
  }))
  by_set <- dplyr::bind_rows(lapply(cols, function(p)
    tibble::tibble(test = sub("^pval_", "", p),
                   set = names(tapply(nc[[p]], nc$set, function(x) 0)),
                   rate_p_lt_alpha = as.numeric(tapply(nc[[p]], nc$set,
                                                       function(x) mean(x < alpha, na.rm = TRUE))))))
  list(per_perm = all, summary = summ, by_set = by_set)
}

gs_donor_scores <- function(objs, sets, celltype_col, donor_col = "patient_ID",
                            min_cells = DEG_MIN_CELLS_PER_SAMPLE, min_genes = 5,
                            include_all = TRUE, min_donors = GS_MIN_DONORS_PER_GROUP) {
  purrr::imap_dfr(objs, function(so, tissue) {
    counts <- tryCatch(SeuratObject::LayerData(so, assay = "RNA", layer = "counts"),
                       error = function(e) SeuratObject::LayerData(
                         SeuratObject::JoinLayers(so[["RNA"]]), layer = "counts"))
    donor  <- as.character(so@meta.data[[donor_col]])
    ct     <- as.character(so@meta.data[[celltype_col]])
    labs   <- sort(unique(stats::na.omit(ct)))
    grp    <- c(if (include_all) list(all = rep(TRUE, ncol(so))),
                stats::setNames(lapply(labs, function(l) ct %in% l), labs))
    purrr::imap_dfr(grp, function(keep, lab) {
      tab <- table(donor[keep]); ok <- names(tab)[tab >= min_cells]
      if (length(ok) < min_donors) return(NULL)
      sel <- which(keep & donor %in% ok)
      f   <- factor(donor[sel], levels = ok)
      mm  <- Matrix::sparse.model.matrix(~ 0 + f); colnames(mm) <- ok
      pb  <- as.matrix(counts[, sel, drop = FALSE] %*% mm)
      dge <- edgeR::DGEList(pb)
      dge <- dge[rowSums(edgeR::cpm(dge) >= 1) >= 2, , keep.lib.sizes = FALSE]
      dge <- edgeR::calcNormFactors(dge)
      lc  <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
      z   <- t(scale(t(lc))); z[!is.finite(z)] <- NA
      purrr::imap_dfr(sets, function(g, s) {
        g <- intersect(g, rownames(z))
        if (length(g) < min_genes) return(NULL)
        tibble::tibble(tissue = tissue, celltype = lab, set = s, donor = ok,
                       score = colMeans(z[g, , drop = FALSE], na.rm = TRUE),
                       n_genes = length(g), n_cells = as.integer(tab[ok]))
      })
    })
  })
}

gs_freeze_params <- function() list(
  GS_MIN_SIZE = GS_MIN_SIZE, GS_MAX_SIZE = GS_MAX_SIZE, GS_MIN_SIZE_SET = GS_MIN_SIZE_SET,
  GS_MATCH_N = GS_MATCH_N, GS_MATCH_BINS = GS_MATCH_BINS, GS_CAMERA_COR_SENS = GS_CAMERA_COR_SENS,
  GS_GSEA_PARAM_SENS = GS_GSEA_PARAM_SENS, GS_IUT_COLS = GS_IUT_COLS,
  GS_MIN_DONORS_PER_GROUP = GS_MIN_DONORS_PER_GROUP, TF_MIN_TARGETS = TF_MIN_TARGETS,
  GS_CUSTOM_RETIRED = sort(GS_CUSTOM_RETIRED), GS_MSIG_EXCLUDE = GS_MSIG_EXCLUDE)

.gs_freeze_text <- function(sets, params) {
  s <- vapply(sort(names(sets)), function(n)
    paste0(n, "\t", paste(sort(unique(sets[[n]])), collapse = ",")), "")
  p <- vapply(sort(names(params)), function(n)
    paste0(n, "=", paste(deparse(params[[n]]), collapse = "")), "")
  c("#params", p, "#sets", s)
}

gs_freeze <- function(file, sets = c(unclass(gs_fetch_msigdb()), GS_CUSTOM),
                      params = gs_freeze_params()) {
  writeLines(.gs_freeze_text(sets, params), file)
  h <- unname(tools::md5sum(file))
  writeLines(h, paste0(file, ".md5"))
  h
}

gs_assert_frozen <- function(file, sets = c(unclass(gs_fetch_msigdb()), GS_CUSTOM),
                             params = gs_freeze_params()) {
  tmp <- tempfile(); on.exit(unlink(tmp))
  writeLines(.gs_freeze_text(sets, params), tmp)
  if (!identical(unname(tools::md5sum(tmp)), readLines(paste0(file, ".md5"), n = 1)))
    stop("gene sets / thresholds differ from ", file, " -- label results exploratory", call. = FALSE)
  invisible(TRUE)
}

tf_get_network <- function(cache = file.path(path$objects, "collectri_human.rds"),
                           csv = TF_NETWORK_CSV) {
  std <- function(n, label) {
    n <- as.data.frame(n)
    if (!"mor" %in% names(n) && "weight" %in% names(n)) n$mor <- n$weight
    if (!nrow(n) || !all(c("source", "target", "mor") %in% names(n)))
      stop("empty or malformed network", call. = FALSE)
    n <- n[, c("source", "target", "mor")]
    attr(n, "network") <- label
    n
  }
  if (!is.null(csv) && file.exists(csv)) return(std(utils::read.csv(csv), paste("csv:", csv)))
  if (file.exists(cache)) return(std(readRDS(cache), "CollecTRI (cached)"))

  try_net <- function(label, expr) tryCatch(std(expr, label), error = function(e) {
    message(label, " failed: ", conditionMessage(e)); NULL })
  net <- NULL
  if (requireNamespace("decoupleR", quietly = TRUE))
    net <- try_net("CollecTRI (decoupleR)",
                   decoupleR::get_collectri(organism = "human", split_complexes = FALSE))
  if (is.null(net) && requireNamespace("OmnipathR", quietly = TRUE))
    net <- try_net("CollecTRI (OmnipathR)", {
      x <- as.data.frame(OmnipathR::collectri())
      data.frame(source = x$source_genesymbol, target = x$target_genesymbol,
                 mor = ifelse(x$is_inhibition == 1 & x$is_stimulation == 0, -1, 1))
    })
  if (is.null(net) && requireNamespace("reticulate", quietly = TRUE) &&
      isTRUE(tryCatch(reticulate::py_module_available("decoupler"), error = function(e) FALSE)))
    net <- try_net("CollecTRI (Python decoupler)", {
      dc  <- reticulate::import("decoupler", convert = FALSE)
      df  <- dc$op$collectri(organism = "human")
      col <- function(nm, type) reticulate::py_to_r(df[[nm]]$astype(type)$tolist())
      data.frame(source = unlist(col("source", "str")),
                 target = unlist(col("target", "str")),
                 mor    = unlist(col("weight", "float")))
    })
  if (!is.null(net)) {
    tryCatch(saveRDS(net, cache), error = function(e)
      warning("could not cache TF network: ", conditionMessage(e), call. = FALSE))
    return(net)
  }
  if (requireNamespace("dorothea", quietly = TRUE)) {
    e <- new.env()
    utils::data("dorothea_hs", package = "dorothea", envir = e)
    d <- e$dorothea_hs[e$dorothea_hs$confidence %in% c("A", "B", "C"), ]
    net <- std(data.frame(source = d$tf, target = d$target, mor = d$mor),
               "DoRothEA A-C (offline fallback)")
    warning("CollecTRI unavailable -- using DoRothEA (confidence A-C) instead. ",
            "Results are labelled; re-run when OmniPath is reachable for CollecTRI.",
            call. = FALSE)
    return(net)
  }
  stop("No TF network: OmniPath unreachable and no offline fallback. ",
       "BiocManager::install('dorothea') for an offline network, or set TF_NETWORK_CSV.",
       call. = FALSE)
}

tf_ulm <- function(stats, net, min_targets = TF_MIN_TARGETS) {
  net <- net[net$target %in% names(stats), , drop = FALSE]
  net <- net[!duplicated(net[, c("source", "target")]), , drop = FALSE]
  nt  <- table(net$source); keep <- names(nt)[nt >= min_targets]
  if (!length(keep)) return(NULL)
  net <- net[net$source %in% keep, , drop = FALSE]
  genes <- names(stats); n <- length(genes)
  X <- Matrix::sparseMatrix(i = match(net$target, genes), j = match(net$source, keep),
                            x = net$mor, dims = c(n, length(keep)))
  y   <- stats - mean(stats)
  Sxx <- Matrix::colSums(X^2) - Matrix::colSums(X)^2 / n
  r   <- as.numeric(Matrix::crossprod(X, y)) / sqrt(Sxx * sum(y^2))
  tv  <- r * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.eps))
  tibble::tibble(tf = keep, score = tv, pval = 2 * stats::pt(-abs(tv), df = n - 2),
                 n_targets = as.integer(nt[keep]),
                 inference = "gene-level, naive p: ranking only")
}

deg_family_report <- function(deg_tbl, comparisons, levels_use, what) {
  d <- deg_tbl[as.character(deg_tbl$comparison) %in% comparisons &
                 as.character(deg_tbl$level) %in% levels_use, , drop = FALSE]
  rep <- d %>%
    dplyr::group_by(tissue, level, subset, comparison) %>%
    dplyr::summarise(rows = dplyr::n(),
                     ok = sum(status %in% "ok"),
                     with_stat = sum(status %in% "ok" & is.finite(stat)),
                     status_first = dplyr::first(status), .groups = "drop") %>%
    dplyr::mutate(reason = dplyr::case_when(
      ok == 0        ~ paste0("fit failed: ", status_first),
      with_stat == 0 ~ "no finite `stat` (column empty or not numeric)",
      TRUE           ~ NA_character_))
  for (cmp in setdiff(comparisons, unique(as.character(rep$comparison))))
    warning(what, ": comparison '", cmp, "' not in the DEG table at level(s) ",
            toString(levels_use), ". Present: ",
            toString(unique(as.character(deg_tbl$comparison))), call. = FALSE)
  bad <- rep[!is.na(rep$reason), , drop = FALSE]
  if (nrow(bad))
    message(what, ": ", nrow(bad), " famil(ies) skipped -- ",
            paste(sprintf("%s | %s | %s | %s: %s", bad$tissue, bad$level, bad$subset,
                          bad$comparison, bad$reason), collapse = "; "))
  rep
}

tf_run <- function(deg_tbl, net, comparisons, levels_use = c("all_myeloid", "celltype")) {
  fam_rep <- deg_family_report(deg_tbl, comparisons, levels_use, "tf_run")
  fam <- fam_rep %>% dplyr::filter(is.na(reason)) %>%
    dplyr::select(tissue, level, subset, comparison)
  out <- dplyr::bind_rows(lapply(seq_len(nrow(fam)), function(i) {
    f  <- fam[i, ]
    df <- deg_tbl %>% dplyr::filter(status == "ok", tissue == f$tissue, level == f$level,
                                    subset == f$subset, comparison == f$comparison)
    r <- tf_ulm(gs_rank_vector(df), net)
    if (is.null(r)) {
      message("tf_run: ", paste(unlist(f), collapse = " | "),
              ": no TF with >= ", TF_MIN_TARGETS, " targets among its genes -- skipped")
      return(NULL)
    }
    r$padj <- stats::p.adjust(r$pval, "BH")
    r$network <- attr(net, "network") %||% NA_character_
    dplyr::bind_cols(f[rep(1, nrow(r)), ], r)
  }))
  for (cmp in setdiff(comparisons, unique(out$comparison)))
    warning("tf_run: no TF results for '", cmp, "' -- see the messages above.", call. = FALSE)
  attr(out, "families") <- fam_rep
  out
}

tf_run_samples <- function(dds, ct, net, min_targets = TF_MIN_TARGETS) {
  cd <- as.data.frame(SummarizedExperiment::colData(dds)); d <- DESeq2::design(dds)
  .gs_check_groups(cd, ct)
  p <- .gs_voom_prep(DESeq2::counts(dds), cd, d, group = ct[1])
  z <- t(scale(t(p$v$E)))
  z <- z[rowSums(!is.finite(z)) == 0, , drop = FALSE]
  sc <- tf_ulm_cells(z, net, min_targets)
  cv <- .gs_contrast_vector(cd, d, ct)
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(sc, p$des), matrix(cv, ncol = 1)))
  nt <- table(unique(net[net$target %in% rownames(z), c("source", "target")])$source)
  tibble::tibble(tf = rownames(sc), score_diff = fit$coefficients[, 1], t = fit$t[, 1],
                 pval = fit$p.value[, 1], padj = stats::p.adjust(fit$p.value[, 1], "BH"),
                 n_targets = as.integer(nt[rownames(sc)]),
                 network = attr(net, "network") %||% NA_character_,
                 inference = "donor-level: limma on per-donor ulm scores")
}

clin_resolve <- function(md, candidates, pattern, what) {
  cols <- if (is.data.frame(md)) colnames(md) else as.character(md)
  hit <- intersect(candidates, cols)
  if (length(hit)) return(hit[1])
  g <- grep(pattern, cols, ignore.case = TRUE, value = TRUE)
  g <- setdiff(g, c("diagnosis", "cerebral_inflammation"))
  if (length(g)) { message(what, ": using column '", g[1], "'"); return(g[1]) }
  message(what, ": no column found (pattern '", pattern, "')")
  NA_character_
}

donor_covariates <- function(objs, cols, donor_col = "patient_ID") {
  keep <- c(donor = donor_col, cols[!is.na(cols)])
  d <- dplyr::bind_rows(lapply(objs, function(so) {
    m <- as.data.frame(lapply(unname(keep), function(k)
      if (k %in% colnames(so@meta.data)) so@meta.data[[k]] else NA),
      stringsAsFactors = FALSE)
    names(m) <- names(keep)
    m[] <- lapply(m, function(v) if (is.factor(v)) as.character(v) else v)
    unique(m)
  }))
  conflict <- d %>% dplyr::group_by(donor) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   ~ dplyr::n_distinct(.x[!is.na(.x)]) > 1), .groups = "drop")
  bad <- conflict$donor[rowSums(as.matrix(conflict[, -1, drop = FALSE])) > 0]
  if (length(bad)) warning("conflicting donor-level values for: ", toString(bad),
                           " -- first value kept", call. = FALSE)
  d <- d %>% dplyr::group_by(donor) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ dplyr::first(.x[!is.na(.x)])),
                     .groups = "drop")
  for (nm in setdiff(names(cols), names(d))) d[[nm]] <- NA
  d
}

kable_save <- function(x, name, ..., dir = path$tables) {
  out <- if (inherits(x, "table") && length(dim(x)) == 2) {
    rn <- names(dimnames(x))[1]
    if (is.null(rn) || !nzchar(rn)) rn <- "row"
    tibble::rownames_to_column(as.data.frame.matrix(x), rn)
  } else as.data.frame(x)
  f <- file.path(dir, paste0(name, ".csv"))
  tryCatch({ readr::write_csv(out, f); message("table saved: ", f) },
           error = function(e) warning("could not save ", f, ": ", conditionMessage(e),
                                       call. = FALSE))
  knitr::kable(x, ...)
}

tf_exclusive_targets <- function(net, tfs) {
  sub <- net[net$source %in% tfs, , drop = FALSE]
  n_src <- table(unique(sub[, c("source", "target")])$target)
  shared <- names(n_src)[n_src > 1]
  out <- sub[!sub$target %in% shared, , drop = FALSE]
  out$source <- paste0(out$source, "_excl")
  out
}

tf_overlap <- function(net, tfs) {
  tg <- lapply(stats::setNames(tfs, tfs), function(t) unique(net$target[net$source == t]))
  present <- tfs[lengths(tg[tfs]) > 0]
  if (length(present) < 2) return(tibble::tibble())
  pr <- utils::combn(present, 2, simplify = FALSE)
  dplyr::bind_rows(lapply(pr, function(p) tibble::tibble(
    tf_a = p[1], tf_b = p[2], n_a = length(tg[[p[1]]]), n_b = length(tg[[p[2]]]),
    shared = length(intersect(tg[[p[1]]], tg[[p[2]]])),
    jaccard = length(intersect(tg[[p[1]]], tg[[p[2]]])) / length(union(tg[[p[1]]], tg[[p[2]]])))))
}

tf_ulm_cells <- function(expr, net, min_targets = TF_MIN_TARGETS, chunk = 5000) {
  genes <- rownames(expr); n <- length(genes)
  net <- net[net$target %in% genes, , drop = FALSE]
  net <- net[!duplicated(net[, c("source", "target")]), , drop = FALSE]
  nt  <- table(net$source); keep <- names(nt)[nt >= min_targets]
  if (!length(keep)) stop("no TF with >= ", min_targets, " targets in the matrix", call. = FALSE)
  net <- net[net$source %in% keep, , drop = FALSE]
  X   <- Matrix::sparseMatrix(i = match(net$target, genes), j = match(net$source, keep),
                              x = net$mor, dims = c(n, length(keep)))
  sx  <- Matrix::colSums(X); Sxx <- Matrix::colSums(X^2) - sx^2 / n
  idx <- split(seq_len(ncol(expr)), ceiling(seq_len(ncol(expr)) / chunk))
  out <- do.call(cbind, lapply(idx, function(j) {
    Y   <- expr[, j, drop = FALSE]
    sy  <- Matrix::colSums(Y); Syy <- Matrix::colSums(Y^2) - sy^2 / n
    Sxy <- as.matrix(Matrix::crossprod(X, Y)) - outer(sx, sy) / n
    r   <- Sxy / sqrt(outer(Sxx, Syy))
    r * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.eps))
  }))
  dimnames(out) <- list(keep, colnames(expr))
  out
}

gs_theme <- function() ggplot2::theme_classic(base_size = FIG_BASE_PT,
                                              base_family = fig_font(quiet = TRUE)) +
  ggplot2::theme(axis.text = ggplot2::element_text(size = FIG_LEGEND_PT, colour = "black"),
                 strip.background = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(size = FIG_LEGEND_PT, face = "bold"),
                 legend.title = ggplot2::element_text(size = FIG_LEGEND_PT),
                 legend.text = ggplot2::element_text(size = FIG_MIN_PT),
                 axis.line = ggplot2::element_line(linewidth = FIG_LINE_MM),
                 axis.ticks = ggplot2::element_line(linewidth = FIG_LINE_MM))

gs_plot_dots <- function(res, provenance, title = NULL, padj_max = 0.05) {
  d <- res %>%
    dplyr::select(-dplyr::any_of(c("panel", "panel.x", "panel.y"))) %>%
    dplyr::left_join(provenance[, c("set", "panel")], by = "set") %>%
    dplyr::mutate(both = if ("agree" %in% names(.)) agree %in% TRUE else
      if ("padj_iut" %in% names(.)) !is.na(padj_iut) & padj_iut < padj_max else
      !is.na(padj_fgsea) & padj_fgsea < padj_max &
        !is.na(padj_camera) & padj_camera < padj_max,
      padj_call = if ("padj_camera_vc_est" %in% names(.)) dplyr::coalesce(padj_camera_vc_est, padj_fgsea) else padj_fgsea,
      set = factor(set, levels = rev(unique(provenance$set[order(provenance$panel)]))),
      subset = factor(subset, levels = unique(c("all", sort(unique(subset))))))
  ggplot2::ggplot(d, ggplot2::aes(subset, set)) +
    ggplot2::geom_point(ggplot2::aes(fill = NES, size = -log10(pmax(padj_call, 1e-10))),
                        shape = 21, colour = "grey80", stroke = 0.2) +
    ggplot2::geom_point(data = d[d$both, ], ggplot2::aes(size = -log10(pmax(padj_call, 1e-10))),
                        shape = 21, fill = NA, colour = "black", stroke = 0.5) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, name = "NES") +
    ggplot2::scale_size_area(max_size = 3, name = "-log10 padj") +
    ggplot2::facet_grid(panel ~ tissue + comparison, scales = "free", space = "free") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    gs_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                                       size = FIG_LEGEND_PT),
                   strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
}

tf_plot_heatmap <- function(res, tfs, title = NULL, padj_max = 0.05) {
  d <- res %>% dplyr::filter(tf %in% tfs) %>%
    dplyr::mutate(tf = factor(tf, levels = rev(tfs)),
                  subset = factor(subset, levels = unique(c("all", sort(unique(subset))))),
                  star = ifelse(padj < padj_max, "*", ""))
  lim <- max(abs(d$score), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(subset, tf, fill = score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = star), size = FIG_LEGEND_PT / .pt, vjust = 0.75) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-lim, lim), name = "Activity\n(ulm t)") +
    ggplot2::facet_grid(. ~ tissue + comparison, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    gs_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                                       size = FIG_LEGEND_PT),
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
}

gs_raw <- function(res) {
  if (!"ranking" %in% names(res)) return(res)
  r <- res[res$ranking == "raw", , drop = FALSE]; r$ranking <- NULL; r
}

GS_TIER_Q_EXPLORATORY <- 0.25
gs_tier <- function(padj_est, p_est, dir_est, padj_fry, dir_fry, legacy = FALSE, padj_max = 0.05) {
  donor <- !is.na(padj_est)
  dplyr::case_when(
    donor & padj_est < padj_max                                       ~ "significant",
    donor & !is.na(padj_fry) & padj_fry < padj_max & p_est < padj_max &
      dir_fry == dir_est                                              ~ "supported",
    donor & padj_est < GS_TIER_Q_EXPLORATORY                          ~ "exploratory",
    !donor & legacy %in% TRUE                                         ~ "screen",
    TRUE                                                              ~ NA_character_)
}
gs_tier_call <- function(tier) tier %in% c("significant", "supported")

gs_finalise <- function(res_all, provenance, res_voom = NULL, padj_max = 0.05) {
  raw <- gs_raw(res_all)
  if (!nrow(raw)) return(raw)
  if ("ranking" %in% names(res_all) && any(res_all$ranking != "raw")) {
    sens <- res_all[res_all$ranking != "raw", c(GS_KEYS, "set", "ranking", "NES", "padj_camera")]
    sens <- tidyr::pivot_wider(sens, names_from = ranking, values_from = c(NES, padj_camera))
    raw  <- dplyr::left_join(raw, sens, by = c(GS_KEYS, "set"))
  }
  ctrl <- raw %>% dplyr::filter(is_control %in% TRUE) %>%
    dplyr::mutate(moved = !is.na(padj_fgsea) & padj_fgsea < padj_max &
                    !is.na(padj_camera) & padj_camera < padj_max) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(GS_KEYS))) %>%
    dplyr::summarise(controls_moved = paste(set[moved], collapse = ","), .groups = "drop") %>%
    dplyr::mutate(global_shift_flag = grepl("RIBOSOME_RP|OXPHOS_SUBUNITS|MITO_RIBOSOME", controls_moved),
                  ieg_flag = grepl("VAN_DEN_BRINK_DISSOCIATION|IEG_DISSOCIATION_STRESS", controls_moved))
  if (!any(raw$set %in% c("VAN_DEN_BRINK_DISSOCIATION", "IEG_DISSOCIATION_STRESS"))) ctrl$ieg_flag <- NA
  raw <- raw %>% dplyr::left_join(ctrl, by = GS_KEYS) %>%
    dplyr::select(-dplyr::any_of("panel")) %>%
    dplyr::left_join(provenance[, c("set", "panel")], by = "set")
  legacy <- !is.na(raw$padj_fgsea) & raw$padj_fgsea < padj_max &
    !is.na(raw$padj_camera) & raw$padj_camera < padj_max &
    sign(raw$NES) == ifelse(raw$camera_direction == "Up", 1, -1)
  if (!is.null(res_voom) && nrow(res_voom)) {
    comb <- gs_combine(res_all, res_voom) %>%
      dplyr::select(dplyr::all_of(GS_KEYS), set, dplyr::any_of(c(
        "mean_logFC", "mean_t", "direction_camera_vc", "pval_camera_vc", "padj_camera_vc",
        "corr_camera_vc_est", "direction_camera_vc_est", "pval_camera_vc_est", "padj_camera_vc_est",
        "direction_fry", "pval_fry", "padj_fry", "p_iut", "padj_iut")))
    raw <- dplyr::left_join(raw, comb, by = c(GS_KEYS, "set"))
  }
  for (cc in c("padj_iut", "padj_camera_vc_est", "direction_camera_vc_est", "padj_fry", "direction_fry"))
    if (!cc %in% names(raw)) raw[[cc]] <- if (grepl("^padj", cc)) NA_real_ else NA_character_
  donor <- !is.na(raw$padj_camera_vc_est)
  raw$inference <- ifelse(donor, "donor-level (voom camera, estimated correlation)", "gene-level screen only")
  raw$tier  <- gs_tier(raw$padj_camera_vc_est, raw$pval_camera_vc_est, raw$direction_camera_vc_est,
                       raw$padj_fry, raw$direction_fry, legacy, padj_max)
  raw$agree <- gs_tier_call(raw$tier)
  raw$fry_sig <- ifelse(donor, !is.na(raw$padj_fry) & raw$padj_fry < padj_max &
                          raw$direction_fry == raw$direction_camera_vc_est, NA)
  raw$robust_matched <- raw$agree & !is.na(raw$padj_matched) & raw$padj_matched < padj_max
  if ("padj_camera_intensity_resid" %in% names(raw))
    raw$robust_intensity <- raw$agree & !is.na(raw$padj_camera_intensity_resid) &
      raw$padj_camera_intensity_resid < padj_max & sign(raw$NES_intensity_resid) == sign(raw$NES)
  raw
}

gs_loo_voom <- function(registry, sets, gs_res, comparisons = DEG_LOO_GS_CONTRASTS,
                        tissues = DEG_LOO_TISSUES, padj_max = 0.05, min_frac_retained = DEG_LOO_STABLE) {
  env <- environment(.gs_check_groups); old <- get("GS_MIN_DONORS_PER_GROUP", envir = env)
  assign("GS_MIN_DONORS_PER_GROUP", max(2L, old - 1L), envir = env)
  on.exit(assign("GS_MIN_DONORS_PER_GROUP", old, envir = env), add = TRUE)
  out <- list()
  for (e in registry) for (cn in names(e$contrasts)) {
    if (!is.null(comparisons) && !cn %in% comparisons) next
    if (!e$tissue %in% tissues) next
    full <- gs_res %>%
      dplyr::filter(tissue == e$tissue, level == e$level, subset == e$subset, comparison == cn,
                    agree %in% TRUE, !is.na(padj_camera_vc_est))
    if (!nrow(full)) next
    ct  <- e$contrasts[[cn]]
    cd  <- as.data.frame(SummarizedExperiment::colData(e$dds))
    cnt <- DESeq2::counts(e$dds); des <- DESeq2::design(e$dds)
    fold <- function(j) tryCatch(.gs_test_voom(cnt[, -j, drop = FALSE], droplevels(cd[-j, , drop = FALSE]),
                                               des, ct, sets),
                                 error = function(err) conditionMessage(err))
    folds <- if (DEG_LOO_CORES > 1 && .Platform$OS.type == "unix")
      parallel::mclapply(seq_len(ncol(cnt)), fold, mc.cores = DEG_LOO_CORES) else lapply(seq_len(ncol(cnt)), fold)
    donors <- sub("__.*$", "", colnames(cnt))
    ok <- vapply(folds, is.data.frame, logical(1))
    n_fail <- sum(!ok); folds <- folds[ok]; names(folds) <- donors[ok]
    if (!length(folds)) next
    pick <- function(r, col) r[[col]][match(full$set, r$set)]
    keep <- sapply(folds, function(r) {
      tr <- gs_tier(pick(r, "padj_camera_vc_est"), pick(r, "pval_camera_vc_est"), pick(r, "direction_camera_vc_est"),
                    pick(r, "padj_fry"), pick(r, "direction_fry"), padj_max = padj_max)
      gs_tier_call(tr) & pick(r, "direction_camera_vc_est") == full$direction_camera_vc_est })
    eff  <- sapply(folds, function(r) pick(r, "mean_logFC"))
    keep <- matrix(keep, nrow = nrow(full), dimnames = list(full$set, names(folds)))
    eff  <- matrix(eff,  nrow = nrow(full), dimnames = list(full$set, names(folds)))
    ratio <- eff / full$mean_logFC
    ps <- tibble::tibble(
      tissue = e$tissue, level = e$level, subset = e$subset, comparison = cn, set = full$set,
      loo_frac_retained    = rowMeans(keep),
      loo_sign_flip        = rowSums(sign(eff) != sign(full$mean_logFC), na.rm = TRUE) > 0,
      loo_min_effect_ratio = apply(ratio, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)),
      loo_lost_without     = apply(!keep, 1, function(z) paste(colnames(keep)[z], collapse = ";")),
      loo_n_folds          = ncol(keep))
    ps$loo_status <- dplyr::case_when(
      n_fail > 0                                ~ sprintf("incomplete (%d/%d folds)", ncol(keep), ncol(keep) + n_fail),
      ps$loo_sign_flip                          ~ "sign_flip",
      ps$loo_frac_retained >= min_frac_retained ~ "stable",
      TRUE                                      ~ "fragile")
    out[[length(out) + 1]] <- ps
    message(sprintf("gene-set LOO | %s | %s | %s | %s: %d sets, %d folds", e$tissue, e$level, e$subset, cn,
                    nrow(full), ncol(keep)))
  }
  dplyr::bind_rows(out)
}

gs_registry_from_units <- function(units,
                                   runner = function(so) run_pseudobulk(so, return_dds = TRUE),
                                   contrasts = list(xALD_vs_ctrl = c("diagnosis", "X-ALD", "control"))) {
  reg <- list()
  for (u in units) {
    so  <- tryCatch(u$get(), error = function(e) NULL)
    dds <- if (is.null(so)) NULL else tryCatch(runner(so), error = function(e) {
      message("gs registry ", u$tissue, " | ", u$level, " | ", u$subset, ": ", conditionMessage(e)); NULL })
    if (!is.null(dds)) reg <- gs_registry_add(reg, u$tissue, u$level, u$subset, dds, contrasts,
                                              qc = gs_sample_qc(so@meta.data))
  }
  reg
}

GS_ADJUST_CDR <- TRUE
GS_CDR_COL    <- "nFeature_RNA"

gs_sample_qc <- function(md, qc_col = GS_CDR_COL, donor_col = "patient_ID", diagnosis_col = "diagnosis") {
  if (!all(c(qc_col, donor_col, diagnosis_col) %in% names(md))) return(NULL)
  key <- paste(as.character(md[[donor_col]]), as.character(md[[diagnosis_col]]), sep = "__")
  tapply(as.numeric(md[[qc_col]]), key, mean, na.rm = TRUE)
}

gs_registry_adjust <- function(reg, name = "cdr_z") {
  out <- lapply(reg, function(e) {
    lab <- paste(e$tissue, e$level, e$subset, sep = " | ")
    q <- if (!is.null(e$qc)) e$qc[colnames(e$dds)] else NULL
    if (is.null(q) || anyNA(q) || any(q <= 0) || stats::sd(log10(q)) == 0) {
      message("gs_registry_adjust: no usable QC for ", lab, " -- skipped"); return(NULL) }
    dds <- e$dds
    dds[[name]] <- as.numeric(scale(log10(q)))
    DESeq2::design(dds) <- stats::update(DESeq2::design(dds), stats::as.formula(paste("~ . +", name)))
    mm <- stats::model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
    if (nrow(mm) - ncol(mm) < 2) {
      message("gs_registry_adjust: too few residual df for ", lab, " -- skipped"); return(NULL) }
    dds <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE), error = function(err) {
      message("gs_registry_adjust: ", lab, ": ", conditionMessage(err)); NULL })
    if (is.null(dds)) return(NULL)
    e$dds <- dds; e$adjusted <- name
    e
  })
  Filter(Negate(is.null), out)
}

gs_registry_deg <- function(reg) {
  dplyr::bind_rows(lapply(reg, function(e) dplyr::bind_rows(lapply(names(e$contrasts), function(cn) {
    r <- tryCatch(DESeq2::results(e$dds, contrast = e$contrasts[[cn]]), error = function(err) NULL)
    if (is.null(r)) return(NULL)
    tibble::tibble(tissue = e$tissue, level = e$level, subset = e$subset, comparison = cn,
                   gene = rownames(r), baseMean = r$baseMean, log2FC_MLE = r$log2FoldChange,
                   stat = r$stat, pvalue = r$pvalue, padj = r$padj, status = "ok")
  }))))
}

gs_deg_concordance <- function(deg_unadj, deg_adj, padj_max = DEG_PADJ_THRESH) {
  dplyr::inner_join(deg_unadj, deg_adj, by = c(GS_KEYS, "gene"), suffix = c("_unadj", "_adj")) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(GS_KEYS))) %>%
    dplyr::summarise(
      n_genes      = dplyr::n(),
      rho_stat     = suppressWarnings(stats::cor(stat_unadj, stat_adj, method = "spearman", use = "complete.obs")),
      n_sig_unadj  = sum(padj_unadj < padj_max, na.rm = TRUE),
      n_sig_adj    = sum(padj_adj < padj_max, na.rm = TRUE),
      n_kept       = sum(padj_unadj < padj_max & padj_adj < padj_max & sign(stat_unadj) == sign(stat_adj), na.rm = TRUE),
      n_sign_flip  = sum(padj_unadj < padj_max & sign(stat_unadj) != sign(stat_adj), na.rm = TRUE),
      n_new        = sum(!(padj_unadj < padj_max) %in% TRUE & padj_adj < padj_max, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::mutate(share_kept = ifelse(n_sig_unadj > 0, n_kept / n_sig_unadj, NA_real_))
}

gs_corr_summary <- function(voom_unadj, voom_adj, padj_max = 0.05) {
  q <- function(x, p) stats::quantile(x, p, na.rm = TRUE, names = FALSE)
  one <- function(v, tag) v %>% dplyr::filter(!is_control) %>%
    dplyr::group_by(tissue, level, comparison) %>%
    dplyr::summarise(n_sets = dplyr::n(), corr_q25 = q(corr_camera_vc_est, .25),
                     corr_median = q(corr_camera_vc_est, .5), corr_q75 = q(corr_camera_vc_est, .75),
                     n_sig_camera_0.01 = sum(padj_camera_vc < padj_max, na.rm = TRUE),
                     n_sig_camera_est  = sum(padj_camera_vc_est < padj_max, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(model = tag)
  dplyr::bind_rows(one(voom_unadj, "unadjusted"), one(voom_adj, "+ detection rate")) %>%
    dplyr::arrange(tissue, level, comparison, model)
}

gs_calibrate_voom <- function(e, sets, ct = e$contrasts[[1]], n_max = 200, alpha = 0.05, seed = 1,
                              control_sets = gs_control_names()) {
  cd  <- as.data.frame(SummarizedExperiment::colData(e$dds))
  cnt <- DESeq2::counts(e$dds); des <- DESeq2::design(e$dds)
  g   <- cd[[ct[1]]]; n <- length(g); k <- sum(g == ct[2]); true_idx <- which(g == ct[2])
  combos <- if (choose(n, k) <= n_max) utils::combn(n, k, simplify = FALSE) else {
    set.seed(seed); unique(replicate(n_max, sort(sample.int(n, k)), simplify = FALSE)) }
  combos <- Filter(function(i) !identical(as.integer(i), as.integer(true_idx)), combos)
  expect <- k * k / n
  d <- dplyr::bind_rows(lapply(seq_along(combos), function(j) {
    lab <- rep(ct[3], n); lab[combos[[j]]] <- ct[2]
    cd2 <- cd; cd2[[ct[1]]] <- factor(lab, levels = levels(g))
    r <- tryCatch(.gs_voom_core(cnt, cd2, des, ct, sets), error = function(err) NULL)
    if (is.null(r)) return(NULL)
    r <- r[!r$set %in% control_sets, , drop = FALSE]
    tibble::tibble(perm = j, overlap = length(intersect(combos[[j]], true_idx)), set = r$set,
                   p_camera_0.01 = r$pval_camera_vc, p_camera_est = r$pval_camera_vc_est, p_fry = r$pval_fry)
  }))
  if (!nrow(d)) return(NULL)
  d$balanced <- abs(d$overlap - expect) < 1
  summ <- dplyr::bind_rows(
    d %>% dplyr::mutate(relabelings = "all"),
    d %>% dplyr::filter(balanced) %>% dplyr::mutate(relabelings = "balanced")) %>%
    dplyr::group_by(relabelings) %>%
    dplyr::summarise(n_relabelings = dplyr::n_distinct(perm), n_tests = dplyr::n(),
                     dplyr::across(dplyr::starts_with("p_"), ~ mean(.x < alpha, na.rm = TRUE),
                                   .names = "rate_{sub('^p_', '', .col)}"),
                     .groups = "drop") %>%
    dplyr::mutate(tissue = e$tissue, level = e$level, subset = e$subset,
                  model = if (is.null(e$adjusted)) "unadjusted" else "+ detection rate",
                  alpha = alpha, .before = 1)
  list(summary = summ, values = d)
}

GS_SET_RENAMED <- c(
  LAM_CORE                   = "KEREN_SHAUL_DAM",
  LYSOSOMAL_DEGRADATIVE      = "KEGG_LYSOSOME",
  LXR_TARGETS                = "REACTOME_NR1H3_NR1H2_REGULATE_GENE_EXPRESSION_LINKED_TO_CHOLESTEROL_TRANSPORT_AND_EFFLUX",
  PEROXISOMAL_BETA_OXIDATION = "REACTOME_BETA_OXIDATION_OF_VERY_LONG_CHAIN_FATTY_ACIDS",
  IEG_DISSOCIATION_STRESS    = "VAN_DEN_BRINK_DISSOCIATION",
  SREBP2_TARGETS_CORE        = "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP",
  MEVALONATE_PRE_SQUALENE    = "KEGG_TERPENOID_BACKBONE_BIOSYNTHESIS",
  STEROL_POST_SQUALENE       = "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  SREBF1_LIPOGENIC           = "REACTOME_NR1H2_NR1H3_REGULATE_GENE_EXPRESSION_LINKED_TO_LIPOGENESIS",
  LXR_TARGETS_STRICT         = "REACTOME_NR1H3_NR1H2_REGULATE_GENE_EXPRESSION_LINKED_TO_CHOLESTEROL_TRANSPORT_AND_EFFLUX",
  INFLAMMASOME_PRIMING       = "REACTOME_THE_NLRP3_INFLAMMASOME",
  NLRP3_INFLAMMASOME         = "REACTOME_THE_NLRP3_INFLAMMASOME",
  VLCFA_PEROXISOME           = "REACTOME_BETA_OXIDATION_OF_VERY_LONG_CHAIN_FATTY_ACIDS",
  ER_STRESS_UPR_CORE         = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  ER_STRESS_UPR_TARGETS      = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
GS_SET_DROPPED <- c("OXYSTEROL_METABOLISM", "LIPID_STORAGE_ESTERIFICATION", "FERROPTOSIS",
                    "VLCFA_ELONGATION", "PHAGOCYTOSIS_EFFEROCYTOSIS", "EFFEROCYTOSIS_RECEPTORS",
                    "MHCII_ANTIGEN_PRESENTATION", "GC_RESPONSE_CANONICAL", "GC_RESPONSE_MYELOID")
gs_current_names <- function(x, quiet = FALSE) {
  y <- ifelse(x %in% names(GS_SET_RENAMED), unname(GS_SET_RENAMED[x]), x)
  drop <- y %in% GS_SET_DROPPED
  if (!quiet) {
    ch <- x != y & !drop
    if (any(ch)) message("gene sets renamed: ", paste(x[ch], "->", y[ch], collapse = "; "))
    if (any(drop)) message("gene sets no longer in the library (use single genes): ", toString(x[drop]))
  }
  unique(y[!drop])
}

gs_shift_diagnostics <- function(deg_tbl, comparisons, levels_use, res_all = NULL) {
  d <- deg_tbl[deg_tbl$status %in% "ok" & as.character(deg_tbl$comparison) %in% comparisons &
                 as.character(deg_tbl$level) %in% levels_use & is.finite(deg_tbl$stat) &
                 is.finite(deg_tbl$baseMean), , drop = FALSE]
  fam <- d %>% dplyr::group_by(dplyr::across(dplyr::all_of(GS_KEYS))) %>%
    dplyr::group_modify(function(x, k) {
      lx <- log10(x$baseMean + 1)
      q  <- cut(rank(lx, ties.method = "first"), 5, labels = FALSE)
      tr <- tryCatch({
        fit <- stats::loess(y ~ lx, data = data.frame(y = x$stat, lx = lx), span = 0.3,
                            family = "symmetric", degree = 1)
        p <- stats::predict(fit, data.frame(lx = stats::quantile(lx, c(0.05, 0.95))))
        unname(p[2] - p[1])
      }, error = function(e) NA_real_)
      tibble::tibble(n_genes = nrow(x), median_stat = stats::median(x$stat),
                     frac_stat_pos = mean(x$stat > 0),
                     median_lfc = stats::median(x$log2FC_MLE, na.rm = TRUE),
                     rho_stat_expr = suppressWarnings(stats::cor(x$stat, lx, method = "spearman")),
                     stat_top_minus_bottom_quintile = stats::median(x$stat[q == 5]) - stats::median(x$stat[q == 1]),
                     trend_p05_p95 = tr)
    }) %>% dplyr::ungroup()
  if (!is.null(res_all) && nrow(res_all)) {
    ra <- if ("ranking" %in% names(res_all)) res_all else dplyr::mutate(res_all, ranking = "raw")
    s <- ra %>% dplyr::filter(!is_control %in% TRUE) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(GS_KEYS)), ranking) %>%
      dplyr::summarise(frac_NES_pos = mean(NES > 0, na.rm = TRUE), mean_NES = mean(NES, na.rm = TRUE),
                       .groups = "drop") %>%
      tidyr::pivot_wider(names_from = ranking, values_from = c(frac_NES_pos, mean_NES))
    ctl <- ra %>% dplyr::filter(ranking == "raw", set %in% c("RIBOSOME_RP", "OXPHOS_SUBUNITS", "MITO_RIBOSOME")) %>%
      dplyr::select(dplyr::all_of(GS_KEYS), set, NES) %>%
      tidyr::pivot_wider(names_from = set, values_from = NES, names_prefix = "NES_")
    fam <- fam %>% dplyr::left_join(s, by = GS_KEYS) %>% dplyr::left_join(ctl, by = GS_KEYS)
  }
  one_sided <- function(f) !is.na(f) & (f <= 0.15 | f >= 0.85)
  balanced  <- function(f) !is.na(f) & f > 0.25 & f < 0.75
  fr  <- if ("frac_NES_pos_raw" %in% names(fam)) fam$frac_NES_pos_raw else NA_real_
  frr <- if ("frac_NES_pos_intensity_resid" %in% names(fam)) fam$frac_NES_pos_intensity_resid else NA_real_
  fam$verdict <- dplyr::case_when(
    one_sided(fr) & balanced(frr) ~ "one-sided sets, balanced after residualising on expression: expression-dependent bias",
    one_sided(fr)                 ~ "one-sided sets, NOT removed by residualising: sets differ from the background as a class",
    abs(fam$rho_stat_expr) >= 0.1 ~ "Wald stat depends on expression level; sets balanced",
    is.na(fr)                     ~ "gene-level only (run gs_run() for the set part)",
    TRUE                          ~ "no clear shift")
  fam
}

gs_plot_shift_ma <- function(deg_tbl, families, max_points = 3000, seed = 1) {
  d <- deg_tbl %>% dplyr::filter(status %in% "ok", is.finite(stat), is.finite(baseMean)) %>%
    dplyr::semi_join(families, by = GS_KEYS) %>%
    dplyr::mutate(family = paste(tissue, subset, comparison, sep = " | "), lx = log10(baseMean + 1))
  pts <- with_seed(seed, d %>% dplyr::group_by(family) %>% dplyr::slice_sample(n = max_points) %>%
                     dplyr::ungroup())
  ggplot2::ggplot(pts, ggplot2::aes(lx, stat)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(size = 0.2, alpha = 0.25, colour = "grey40") +
    ggplot2::geom_smooth(data = d, method = "loess", formula = y ~ x, span = 0.3, se = FALSE,
                         colour = "#B2182B", linewidth = 0.5) +
    ggplot2::facet_wrap(~ family, scales = "free_y") +
    ggplot2::labs(x = "log10 (baseMean + 1)", y = "Wald statistic",
                  title = "Wald statistic along expression level (red: loess; flat at 0 = no expression-dependent shift)") +
    gs_theme()
}

gs_donor_shift_qc <- function(objs, celltype_col, donor_col = "patient_ID", group_col = "diagnosis",
                              qc_cols = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent_mt",
                                          "pct_mt", "percent.ribo", "percent_ribo"),
                              include_celltypes = FALSE, min_cells = DEG_MIN_CELLS_PER_SAMPLE) {
  universe <- unique(unlist(lapply(objs, rownames)))
  sc <- gs_donor_scores(objs, gs_control_sets(universe), celltype_col, donor_col,
                        min_cells = min_cells, include_all = TRUE)
  if (!include_celltypes) sc <- sc[sc$celltype == "all", , drop = FALSE]
  sc$n_cells <- NULL
  qc <- purrr::imap_dfr(objs, function(so, ti) {
    md <- so@meta.data
    q  <- intersect(qc_cols, names(md))
    ct <- as.character(md[[celltype_col]])
    labs <- sort(unique(stats::na.omit(ct)))
    grp <- c(list(all = rep(TRUE, nrow(md))),
             if (include_celltypes) stats::setNames(lapply(labs, function(l) ct %in% l), labs))
    purrr::imap_dfr(grp, function(keep, lab) {
      md[keep, , drop = FALSE] %>%
        dplyr::group_by(donor = as.character(.data[[donor_col]])) %>%
        dplyr::summarise(group = dplyr::first(as.character(.data[[group_col]])), n_cells = dplyr::n(),
                         dplyr::across(dplyr::all_of(q), ~ stats::median(as.numeric(.x), na.rm = TRUE),
                                       .names = "median_{.col}"), .groups = "drop") %>%
        dplyr::mutate(tissue = ti, celltype = lab)
    })
  })
  d <- dplyr::inner_join(sc, qc, by = c("tissue", "celltype", "donor"))
  qcm <- c("n_cells", grep("^median_", names(d), value = TRUE))
  wmean <- function(z, g) stats::ave(z, g, FUN = function(v) mean(v, na.rm = TRUE))
  summ <- d %>% tidyr::pivot_longer(dplyr::all_of(qcm), names_to = "qc", values_to = "qc_value") %>%
    dplyr::group_by(tissue, celltype, set, qc) %>%
    dplyr::summarise(
      n_donors = dplyr::n(),
      score_diff_xald_ctrl = mean(score[group != "control"], na.rm = TRUE) - mean(score[group == "control"], na.rm = TRUE),
      qc_diff_xald_ctrl_sd = (mean(qc_value[group != "control"], na.rm = TRUE) -
                                mean(qc_value[group == "control"], na.rm = TRUE)) / stats::sd(qc_value, na.rm = TRUE),
      rho_all = suppressWarnings(stats::cor(score, qc_value, method = "spearman", use = "complete.obs")),
      rho_within_group = suppressWarnings(stats::cor(score - wmean(score, group), qc_value - wmean(qc_value, group),
                                                     method = "spearman", use = "complete.obs")),
      slope_within = {
        y <- score - wmean(score, group); x <- qc_value - wmean(qc_value, group)
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) >= 3 && stats::var(x[ok]) > 0) unname(stats::coef(stats::lm(y[ok] ~ x[ok]))[2]) else NA_real_
      },
      qc_diff_raw = mean(qc_value[group != "control"], na.rm = TRUE) - mean(qc_value[group == "control"], na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::mutate(
      share_explained = slope_within * qc_diff_raw / score_diff_xald_ctrl,
      reading = dplyr::case_when(
        n_donors < 6 ~ "",
        abs(rho_within_group) >= 0.6 & abs(qc_diff_xald_ctrl_sd) >= 0.5 & share_explained >= 0.5 ~
          "follows QC within groups and QC differs between groups: may explain the shift (technical candidate)",
        abs(rho_within_group) >= 0.6 & abs(qc_diff_xald_ctrl_sd) < 0.5 ~
          "follows QC within groups, but QC barely differs between groups: adds noise, does not explain the shift",
        abs(rho_within_group) >= 0.6 ~
          "follows QC within groups; the QC difference between groups explains < 50% of the shift or runs against it",
        abs(qc_diff_xald_ctrl_sd) >= 0.5 ~
          "QC differs between groups but does not track the score: not the explanation",
        TRUE ~ "")) %>%
    dplyr::select(-qc_diff_raw)
  list(donor_values = d, summary = summ)
}
