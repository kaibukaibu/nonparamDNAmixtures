###### Functions to use for the simulations


# For turning wide profile data to a long format
wide_to_long <- function(trace){
  n_siballeles <- (ncol(trace)-2)/2
  trace <- trace %>% 
    as.matrix() %>% 
    cbind(matrix(nrow=nrow(trace),
                 ncol=n_siballeles))
  
  for(y in 1:n_siballeles){
    trace[,2+2*n_siballeles+y] <- paste0(trace[,2+y], "__", as.character(as.numeric(trace[,2+n_siballeles+y])))
  }
  
  trace %>% 
    as.data.frame %>% 
    select(!starts_with(c("Allele", "Height"))) %>%
    data.table::as.data.table() %>%
    data.table::melt(id.vars=1:2) %>%
    filter(!str_detect(value, "NA")) %>%
    group_by(SampleName, Marker, value) %>%
    mutate(Allele=str_split(value, "__")[[1]][1],
           Height=str_split(value, "__")[[1]][2]) %>%
    ungroup() %>%
    select(-variable, -value)
  
}





# Version 2
wide_to_long2 <- function(trace){
  trace %>% 
    select(-any_of(starts_with("Size"))) %>%
    select(-any_of(starts_with("Dye"))) %>%
    filter(!(Marker %in% c("AMEL", "Yindel", "DYS391"))) %>% 
    mutate(across(starts_with(c("Allele", "Height")), as.character)) %>%
    pivot_longer(
      cols = starts_with("Allele") | starts_with("Height"),
      names_to = c(".value", "set"),
      names_pattern = "(Allele|Height)(\\d+)"
    ) %>% 
    filter((!Allele %in% c("OL", "")) & !is.na(Allele)) %>% 
    select(-set) %>% 
    arrange(SampleName, Marker, Allele)
}





# For turning long to wide format
long_to_wide <- function(datalong, samplelabel){
  datalong <- datalong %>% 
    rowwise() %>% 
    mutate(Allele_Height = paste0(Allele, "__", Height)) %>% 
    ungroup()
  
  data03 <- data.frame()
  n_alleles_max <- datalong %>% group_by(Marker) %>% summarise(n=n()) %>% pull(n) %>% max()
  for(marker in unique(datalong$Marker)){
    data02 <- datalong %>% 
      filter(marker == Marker) %>% 
      arrange(as.numeric(Allele)) %>% 
      mutate(Allele_Height = factor(Allele_Height, levels=Allele_Height)) %>% 
      data.table::as.data.table() %>% 
      data.table::dcast(formula = Marker ~ Allele_Height, id=c("Sample.File", "Marker"), value.var = "Allele_Height")
    
    if(ncol(data02)-1 != n_alleles_max){
      data02 <- cbind(data02, matrix(rep(NA, n_alleles_max-(ncol(data02)-1)), nrow = 1))
    }
    
    colnames(data02) <- c("Marker", paste0("combo", seq(1:n_alleles_max)))
    data03 <- rbind(data03, data02)
  }
  
  data03 <- data03 %>% 
    mutate(SampleName = samplelabel) %>% 
    relocate(SampleName)
  
  alleles <- apply(as.matrix(data03[,!1:2]), c(1, 2), function(x) str_split(x, "__")[[1]][1])
  heights <- apply(as.matrix(data03[,!1:2]), c(1, 2), function(x) str_split(x, "__")[[1]][2])
  
  data04 <- cbind(data03[,"SampleName"], data03[,"Marker"], alleles, heights)
  colnames(data04) <- c("SampleName",  "Marker", paste0("Allele", 1:n_alleles_max), paste0("Height", 1:n_alleles_max))
  
  data04 %>% 
    mutate(SampleName = gsub(pattern="\\.", replacement="", x=SampleName)) %>% 
    mutate(SampleName = gsub(pattern="_", replacement="-", x=SampleName)) %>% 
    mutate(across(everything(), ~ ifelse(is.na(.), "", .)))
}





