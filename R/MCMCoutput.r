## Plotting posteriors
test1 <- mcmcOutput %>% filter(V1 < mean(V1)) %>%
						 top_n(800) %>%
						 mutate(country = "BE")
						 
test2 <- mcmcOutput %>% filter(V1 >= mean(V1)) %>%
						 top_n(800) %>%
						 mutate(country = "FR")
 
test <- bind_rows(test1, test2) %>%
		rename(gamma0 = V1) %>%
		select(-V2)
		
# Plot density for comparison
ggplot(test, aes(x = gamma0, y = country)) + geom_halfeyeh(.width = c(0.5, 0.95))
 
 
 ## Plotting prevalence
 
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
			   