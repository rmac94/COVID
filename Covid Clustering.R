library(data.table) # data handling
library(TSrepr) # ts representations, use dev version devtools::install_github("PetoLau/TSrepr")
library(ggplot2) # visualisations
library(dtwclust) # clustering using DTW distance
library(dygraphs) # interactive visualizations
library(ggrepel) # nice labels in ggplot
library(dendextend) # dendrograms
library(DT) # nice datatable
library(htmlwidgets) # save html objects
library(keras) # normalize objectsd
library(smooth)



# Pull population estimates by UA -----------------------------------------


#Using MonstR API for ONS data

#Installing and loading the monstR package
remotes::install_github("HFAnalyticsLab/monstR", build_vignettes = TRUE )
library(monstR)

##Choosing a dataset

#The dataset id is key. Keywords are present for searching.  

#ons_available_datasets function will list all available datasets. Choose the id which seems most suitable. 
all_dataset_metadata = ons_available_datasets()

#You can query the editions and versions available for your preferred dataset
ons_available_editions(id = "mid-year-pop-est")
ons_available_versions(id = "mid-year-pop-est", edition = "mid-2019-april-2020-geography")


#The following pipeline downloads the chosen dataset by id into your working directory

datasets <- ons_datasets_setup(monstr_pipeline_defaults())
latest <- datasets %>% ons_dataset_by_id("mid-year-pop-est", edition="mid-2019-april-2020-geography")

#The following is a hacky workaround to read the csv directly into the R Studio environment. Avoid continued usage of this. 
chosen_dataset = read.csv(latest$downloads[1]$csv$href)

pop_process <- chosen_dataset %>% filter(time == 2019 &
                                                       sex == "All" &
                                                       age == "Total") %>%
  select(mnemonic = admin.geography, pop = v4_0)


#UK population
# pop <- read_csv("population estimates by la - mid-2018 - cleaned.csv")
# 
# pop_process <- pop %>% rename(la = `local authority: county / unitary (as of April 2019)`) %>%
#   gather(key = "age range", value = "pop", c(-sex, -la, -mnemonic)) %>%
#   mutate(`age range` = str_sub(`age range`, 1, 5)) %>%
#   mutate(`age range` = ifelse(`age range` == "80+_1", "80+", `age range`)) %>%
#   group_by(mnemonic) %>%
#   summarise(pop = sum(pop))

#UK data load
#CV_regional_source <- read_csv("https://raw.githubusercontent.com/odileeds/coronavirus-data/master/data/coronavirus-cases.csv") 

CV_regional_source <- read_csv("https://coronavirus.data.gov.uk/downloads/csv/coronavirus-cases_latest.csv")

# Establish cumulative time series and join population data
CV_regional_pop <- CV_regional_source %>%
  group_by(`Area name`,`Area code`) %>% 
  arrange(`Area name`, `Specimen date`) %>%
  filter(!is.na(`Cumulative lab-confirmed cases`)) %>%
  select(`Area Code` = `Area code`, `Area Name`=`Area name`, date = `Specimen date`, Cumulative_Cases = `Cumulative lab-confirmed cases`, Daily_total = `Daily lab-confirmed cases`) %>%
  complete(date = seq.Date(min(date), max(date), by="day")) %>%
  fill(Cumulative_Cases) %>%
  mutate(Daily_total = case_when(is.na(Daily_total)~0, 
                                 T ~ Daily_total)) %>%
  left_join(pop_process, by = c("Area Code"="mnemonic")) %>%
  filter(!is.na(pop)) %>%
  distinct()


# Build normalise function 
normalise <- function(x){
  z <- (x - min(x))/(max(x)-min(x))
  return(z)
}


# Sort order and standardise for population size
CV_regional_process <- CV_regional_pop %>%
  arrange(`Area Name`, date) %>%
  #determine standardised cases & take differences
  mutate(Active_cases = case_when(is.na(lag(Daily_total,14)) ~ Cumulative_Cases - 0,
                                      T ~ Cumulative_Cases - lag(Daily_total, 14)),
         delta_cases = case_when(is.na(lag(Active_cases)) ~ Active_cases - 0,
                                 T ~ Active_cases - lag(Active_cases,1)),
         Active_cases_p_000_000 = Active_cases*1000000/pop,
         delta_cases_p_000_000 = delta_cases*1000000/pop,
         Cumulative_Cases_p_000_000 = Cumulative_Cases*1000000/pop) %>%
  #filter for unilateral trajectory starting point
  filter(Cumulative_Cases_p_000_000 > 49) %>%
  # create dummy trajectory date
  mutate(day_from_filter = 1:n(),
  #Calculate moving average over 3 days
  sma_delta_cases_p_000_000 = case_when(is.na(lag(delta_cases_p_000_000,1)) ~ delta_cases_p_000_000,
                                        is.na(lag(delta_cases_p_000_000,2)) ~ (delta_cases_p_000_000 + lag(delta_cases_p_000_000,1))/2,
                                        T ~ (delta_cases_p_000_000 + lag(delta_cases_p_000_000,1) + lag(delta_cases_p_000_000,2))/3)) %>%
  ungroup() %>%
  mutate(norm_sma_delta = normalise(sma_delta_cases_p_000_000)) %>%
  select(-Cumulative_Cases_p_000_000,
         -Active_cases, 
         -pop,
         -delta_cases,
         -date)



#Establish lags and then sequences

n_steps <- 15
tim_ser <- CV_regional_process

segmenting <- vector("list")

for(i in 1:length(unique(tim_ser$`Area Name`))){
  
  Area_subset <- tim_ser %>% filter(`Area Name` == unique(tim_ser$`Area Name`)[i])

  for(j in 1:(max(Area_subset$day_from_filter)-n_steps)){
    Time_subset <- Area_subset %>%
      filter(day_from_filter >= j & day_from_filter < j + n_steps) %>%
      mutate(lag = j)
    segmenting[[paste0(unique(tim_ser$`Area Name`)[i],j)]] <- Time_subset
  }
}

segmenting_output <- segmenting[[1]]
for(i in 1:(length(segmenting)-1)){

data_1 <- segmenting[[i+1]]

segmenting_output <- rbind(segmenting_output, data_1)
}

# Label 14 day splits
segmenting_processed <- segmenting_output %>%
  ungroup() %>%
  mutate(Split = paste0(segmenting_output$`Area Name`, "-", segmenting_output$lag)) %>%
  group_by(`Area Code`,`Area Name`, Split) %>%
  mutate(day_from_segment = 1:n()) %>%
  ungroup()
         
         
# Plot overall splits
ggplot(segmenting_processed) +
  geom_line(aes(#y = delta_cases_p_000_000,
                y = norm_sma_delta,
                x = day_from_segment,
                group = Split),
            alpha = 0.01,
            show.legend = F)


trajectories_list <- vector("list")

for(i in 1:length(unique(segmenting_processed$Split))){
  ts <- segmenting_processed %>%
    filter(Split == unique(segmenting_processed$Split)[i]) %>%
    select(day_from_segment,norm_sma_delta)
  trajectories_list[[i]] <- ts
  
}
names(trajectories_list) <- unique(segmenting_processed$Split)

#Determine number of clusters
n_clusters = 20

# clustering with dtwclust package function
cluster_obj <- tsclust(trajectories_list,
                  type = "hierarchical",
                  k = n_clusters,
                  distance = "dtw_basic",
                  centroid = dba,
                  trace = FALSE,
                  seed = 54321,
                  control = hierarchical_control(method = "ward.D2"),
                  args = tsclust_args(dist = list(norm = "L2"))
)


#Extact cluster ID
cluster_id <- data.frame(
  Cluster = cluster_obj@cluster,
  Split = names(cluster_obj@cluster)
)

#Identify cluster by Cluster ID
cluster_output <- segmenting_processed %>%
  left_join(cluster_id, by = c("Split"="Split"))

#Plot Clusters
ggplot(cluster_output) +
  geom_line(aes(
    y = Active_cases_p_000_000,
    #y = norm_sma_delta,
    x = day_from_segment,
    group = Split),
    alpha = 0.01,
    show.legend = F) +
  facet_wrap(.~Cluster)



