library(tidyverse)

#*Load Data
individuals_all <- read_rds("~/hla_pop_project/scripts/final_scripts/total_run_FINAL/unique_individuals_all_v2")
countries <- read_rds("~/hla_pop_project/scripts/countries")

#Reformat input data:
colnames(individuals_all) <- c("individual", "Epitopes")

#Replace space with underscore in country names:
countries <- gsub(" ", "_", countries)
fix <- individuals_all[,1]
fixed <- gsub(" ", "_", fix)
individuals_all[,1] <- fixed

#Get Country Set Cover Solutions:-----------------------------------------------------------------------------------------------

step.1 <- function(x, S, I) {
  test <- S[[x]]
  test <- test[!(test %in% I)]
  return(test)
}

step.2 <- function(x, S) {length(S[[x]])}

cover_set <- function(country, n) {
  print(country)
  pop.individuals <- individuals_all[grepl(country, individuals_all[,1]),]
  if(NROW(pop.individuals) > 0) {
    pop.individuals <- as.data.frame(pop.individuals)
    
    #Initialize Everything:
    U <- unique(pop.individuals$individual)
    I <- c()
    ep.order <- c()
    ep.sets <- c()
    ep.added <- c()
    
    #Create S, vector of sets:
    pre.S <- pop.individuals %>%
      group_by(Epitopes) %>%
      summarise(s_i = paste(individual, collapse = " "))
    
    S <- lapply(pre.S$s_i, function(x) {unlist(strsplit(x, split = " "))})
    names(S) <- pre.S$Epitopes
    
    
    #Run Algorithm:
    while (!setequal(I, U)) {
      # PART 1: Find S_i with largest length of S_i - I
      test.1 <- lapply(c(1:length(S)), step.1, S, I)
      test.2 <- sapply(c(1:length(S)), step.2, test.1)
      
      num <- which(test.2==max(test.2))[1]
      name.num <- names(S[num]) #Epitope String
      
      # PART 2: Add this S_i's elements to I
      test.3 <- S[name.num][[1]]
      I <- unique(unlist(c(I, test.3)))
      
      # PART 3: Remove S_i from S
      S <- S[-num]
      
      # PART 4: Save Info
      ep.order <- c(ep.order, name.num)
      ep.sets <- c(ep.sets, test.3)
      ep.added <- c(ep.added, length(I))
      
      print(length(I))
    }
    
    
    #Create DF Output:
    pop.setSolution <- tibble(Epitopes_order = ep.order, Total_people = ep.added)
    pop.setSolution <- pop.setSolution %>% mutate(perc = (Total_people/n)*100)
    
    pop.setSolution <- pop.setSolution[order(pop.setSolution$Total_people),]
    pop.setSolution$Epitopes_order <- factor(pop.setSolution$Epitopes_order, 
                                             levels = unique(pop.setSolution$Epitopes_order))
    
    
    return(pop.setSolution)
  } else {
    return(NA)
  }
}


country_setCovers <- lapply(countries, cover_set, n = 100000)
names(country_setCovers) <- countries

write_rds(country_setCovers, "~/hla_pop_project/scripts/final_scripts/total_run/country_setCovers")

