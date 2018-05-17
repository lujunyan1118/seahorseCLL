# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Function to run enrichment analysis in R
#'
#' A funtion to perform GSEA or PAGE analysis
#' @param inputTab Ranked gene list for enrichment analysis, gene symbols as rownames and column stat contains the statistics for the ranking.
#' @param gmtFile A path to the gene signature file.
#' @param GSAmethod Method for enrichment analysis, currently supports GSEA and PAGE
#' @param nPerm Number of permutations for GSEA analysis
#' @export

runGSEA <- function(inputTab,gmtFile,GSAmethod="gsea",nPerm=1000){
  require(piano)
  inGMT <- loadGSC(gmtFile,type="gmt")
  rankTab <- inputTab[order(inputTab[,1],decreasing = TRUE),,drop=FALSE] #re-rank by score

  if (GSAmethod == "gsea"){
    #readin geneset database
    #GSEA analysis
    res <- runGSA(geneLevelStats = rankTab,geneSetStat = GSAmethod,adjMethod = "fdr", gsc=inGMT, signifMethod = 'geneSampling', nPerm = nPerm)
    GSAsummaryTable(res)
  } else if (GSAmethod == "page"){
    res <- runGSA(geneLevelStats = rankTab,geneSetStat = GSAmethod,adjMethod = "fdr", gsc=inGMT, signifMethod = 'nullDist')
    GSAsummaryTable(res)
  }
}

#' Barplot for enrichment analysis result
#'
#' A function to plot enrichment analysis result from runGSEA function
#' @param resTab A list object of result tables from runGSEA function
#' @param pCut The cutoff for p-values of signatures to be plot
#' @param ifFDR Whether to use FDR cutoff or raw p value cutoff (default)
#' @param setName The y axis label of the bar plot
#' @export
#'
plotEnrichmentBar <- function(resTab, pCut = 0.05, ifFDR = FALSE, setName = "Signatures") {
  pList <- list()
  rowNum <- c()
  for (i in names(resTab)) {
    plotTab <- resTab[[i]]
    if (ifFDR) {
      plotTab <- dplyr::filter(plotTab, `p adj (dist.dir.up)` <= pCut | `p adj (dist.dir.dn)` <= pCut)
    } else {
      plotTab <- dplyr::filter(plotTab, `p (dist.dir.up)` <= pCut | `p (dist.dir.dn)` <= pCut)
    }
    if (nrow(plotTab) == 0) {
      print("No sets passed the criteria")
      next
    } else {
      #firstly, process the result table
      plotTab <- apply(plotTab, 1, function(x) {
        statSign <- as.numeric(x[3])
        data.frame(Name = x[1], p = as.numeric(ifelse(statSign >= 0, x[4], x[6])), geneNum = ifelse(statSign >= 0, x[8], x[9]),
                   Direction = ifelse(statSign > 0, "Up", "Down"), stringsAsFactors = FALSE)
      }) %>% do.call(rbind,.)

      plotTab$Name <- sprintf("%s (%s)",plotTab$Name,plotTab$geneNum)
      plotTab <- plotTab[with(plotTab,order(Direction, p, decreasing=TRUE)),]
      plotTab$Direction <- factor(plotTab$Direction, levels = c("Down","Up"))
      plotTab$Name <- factor(plotTab$Name, levels = plotTab$Name)
      #plot the barplot
      pList[[i]] <- ggplot(data=plotTab, aes(x=Name, y= -log10(p), fill=Direction)) +
        geom_bar(position="dodge",stat="identity", width = 0.5) +
        scale_fill_manual(values=c(Up = "blue", Down = "red")) +
        coord_fixed(ratio = 0.5) + coord_flip() + xlab(setName) +
        ylab(expression(-log[10]*'('*p*')')) +
        ggtitle(i) + theme_bw() + theme(plot.title = element_text(face = "bold", hjust =0.5),
                                        axis.title = element_text(size=15))
      rowNum <-c(rowNum,nrow(plotTab))
    }
  }

  if (length(pList) == 0) {
    print("Nothing to plot")
  } else {
    rowNum <- rowNum
    grobList <- lapply(pList, ggplotGrob)
    grobList <- do.call(gridExtra::gtable_rbind,c(grobList,size="max"))
    panels <- grobList$layout$t[grep("panel", grobList$layout$name)]
    grobList$heights[panels] <- unit(rowNum, "null")
  }
  return(grobList)
}
