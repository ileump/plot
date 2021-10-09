trace(ropls::plot, tracer = browser, signature = 'opls')

ropls::plot(pls, typeVc = 'x-score', parDevNewL = F)

# debugonce(.plotF)
# continue
debugonce(.plotLegendF)
# continue

## row 196 in ropls::plot
# for (ploC in typeVc) .plotF(ploC, opl = x, obsColVc = obsColVc, 
#                             obsLabVc = obsLabVc, obsLegVc = obsLegVc, layL = layL, 
#                             parCexN = parCexN, parEllipsesL = parEllipsesL, parTitleL = parTitleL, 
#                             parCompVi = parCompVi, typeVc = typeVc, tCompMN = tCompMN, 
#                             pCompMN = pCompMN, cxtCompMN = cxtCompMN, cytCompMN = cytCompMN, 
#                             topLoadMN = topLoadMN, pexVi = pexVi, tesColVc = tesColVc, 
#                             tesLabVc = tesLabVc, tesLegVc = tesLegVc)
isBlue <- obsColVc == 'blue'
obsColVc[isBlue] <- 'red'
obsColVc[!isBlue] <- 'blue'

## row 71 in .plotLegendF
# text(xLegN, seq(yBotN + (yTopN - yBotN)/(2 * length(scaVc)), 
#                 yTopN - (yTopN - yBotN)/(2 * length(scaVc)), length = length(scaVc)), 
#      adj = c(xAdjN, 0.5), cex = txtCexN, col = scaVc, 
#      labels = names(scaVc))
names(scaVc) <- rev(names(scaVc))

untrace(ropls::plot, signature = 'opls')
