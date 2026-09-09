# ============================================================
# Proksee/CGView-style circular genome map for LyasA and LyasB, built
# with circlize -- each contig is its own sector (circlize's native
# handling of multi-chromosome/multi-contig genomes), arranged by
# size, with tracks for mobileOG hits (colored by category), CARD
# resistance hits, and CRISPR array positions.
# ============================================================

library(circlize)
library(jsonlite)
library(scales)

data <- fromJSON(path.expand("~/Downloads/circular_map_data.json"), simplifyVector = FALSE)

# mobileOG category palette -- reusing project hues where sensible,
# extended for the 5 distinct categories present in this data
category_colors <- c(
  "replication/recombination/repair" = "#B491E5",  # violet
  "stability/transfer/defense"        = "#E4749B",  # pink
  "phage"                              = "#66A8E5",  # sky blue
  "transfer"                            = "#F0C37D",  # amber
  "integration/excision"                 = "#A6B954"   # olive
)

build_map <- function(genome_name, gdata, out_path) {
  contig_lengths <- gdata$contig_lengths
  # Order contigs by size, descending (Proksee convention)
  contig_names <- names(contig_lengths)[order(-unlist(contig_lengths))]
  contig_lens <- unlist(contig_lengths)[contig_names]

  sectors <- data.frame(
    sector = contig_names,
    start = 0,
    end = contig_lens
  )

  png(out_path, width = 3000, height = 3000, res = 300)
  circos.clear()
  circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.6, start.degree = 90)
  circos.initialize(sectors = sectors$sector, xlim = cbind(sectors$start, sectors$end))

  # Track 1: contig backbone (alternating shade for visual separation)
  circos.track(ylim = c(0, 1), track.height = 0.08, bg.border = NA, panel.fun = function(x, y) {
    sector = get.cell.meta.data("sector.index")
    idx = which(contig_names == sector)
    col = if (idx %% 2 == 0) "grey75" else "grey55"
    circos.rect(0, 0, get.cell.meta.data("xlim")[2], 1, col = col, border = NA)
    # label only the largest few contigs -- too many to label all
    if (idx <= 8) {
      circos.text(get.cell.meta.data("xcenter"), 1.5, sector, cex = 0.5, facing = "clockwise",
                    niceFacing = TRUE, adj = c(0, 0.5))
    }
  })

  # Track 2: mobileOG hits, colored by category
  circos.track(ylim = c(0, 1), track.height = 0.12, bg.border = "grey90", panel.fun = function(x, y) {})
  for (hit in gdata$mobileog_hits) {
    if (!is.null(hit$contig) && hit$contig %in% contig_names) {
      col = category_colors[hit$category]
      if (is.na(col)) col = "grey50"
      circos.rect(hit$start, 0, hit$end, 1, sector.index = hit$contig, track.index = 2,
                    col = col, border = col)
    }
  }

  # Track 3: CARD resistance hits
  circos.track(ylim = c(0, 1), track.height = 0.08, bg.border = "grey90", panel.fun = function(x, y) {})
  for (hit in gdata$card_hits) {
    if (!is.null(hit$contig) && hit$contig %in% contig_names) {
      circos.rect(hit$start, 0, hit$end, 1, sector.index = hit$contig, track.index = 3,
                    col = "red", border = "red")
    }
  }

  # Track 4: CRISPR hits
  circos.track(ylim = c(0, 1), track.height = 0.08, bg.border = "grey90", panel.fun = function(x, y) {})
  for (hit in gdata$crispr_hits) {
    if (!is.null(hit$contig) && hit$contig %in% contig_names) {
      circos.rect(hit$start, 0, hit$end, 1, sector.index = hit$contig, track.index = 4,
                    col = "black", border = "black")
    }
  }

  title(paste0(genome_name, " -- contigs, mobileOG hits, CARD resistance, CRISPR"), cex.main = 1.1)

  legend("bottomleft", legend = names(category_colors), fill = category_colors,
         title = "mobileOG category", cex = 0.55, bty = "n")
  legend("bottomright", legend = c("CARD resistance", "CRISPR array"),
         fill = c("red", "black"), cex = 0.6, bty = "n")

  circos.clear()
  dev.off()
  cat("Saved:", out_path, "\n")
}

for (genome in names(data)) {
  out_path <- path.expand(paste0("~/Downloads/", genome, "_circular_map.png"))
  build_map(genome, data[[genome]], out_path)
}
