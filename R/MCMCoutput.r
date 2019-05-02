 results <- FoI(0.141,otherParams,seroData)$prev

 ggplot() +
    geom_point(data = as.data.frame(cbind(age = as.numeric(row.names(htab)), 
                                          prop = htab[, 2] / rowSums(htab),
                                          tot = rowSums(htab))),
               mapping = aes(x = age, y = prop, size = tot),
               pch = 1) +
	geom_line(data = data.frame(unique(cbind(seroData$AGE, results))),
              mapping = aes(x = unique(seroData$AGE), y = results), colour = 4) +
	xlim(0,50) +
	ylim(0,1)
			   