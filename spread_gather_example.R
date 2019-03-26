# These are all the commands that are present in the Rmarkdown file


# install.packages("tidyverse") # you may need to run this

library(tidyverse)

###
CD3 <- read_csv('CD3pos.csv')

CD3[1:3,1:3]
###


###
CD3$pignum <- sub('Specimen_([0-9]+) ([A-Za-z]+)_.*','\\1', CD3$Sample)  # some regular expressions to extract metadata from sample names

CD3$tissue <- sub('Specimen_([0-9]+) ([A-Za-z]+)_.*','\\2', CD3$Sample)

CD3$treatment <- ifelse(CD3$pignum %in% c(67:73), 'control', 'RPS')      # adding in treatment groups...

CD3 <- CD3 %>% select(pignum, tissue, treatment, everything(), -Sample)  # ordering the dataframe so metadata is first, also dropping "Sample" column

CD3
###



# convert to long format (), This data would be considered tidy?
CD3 %>% gather(key = 'cell_type', value = 'num_events', -c(1:3))

CD3.long <- CD3 %>% gather(key = 'cell_type', value = 'num_events', -c(1:3))


######


CD3.long <- CD3.long %>% group_by(pignum, tissue) %>%
  mutate(percent_CD3 = (num_events / sum(num_events)) * 100) %>%
  select(-num_events) %>% ungroup()


#######

CD3.long <- CD3.long %>% mutate(cell_type=gsub('/','', cell_type)) %>% 
  mutate(cell_type=gsub('CD3\\+','', cell_type))

#########


CD3.long %>% ggplot(aes(x=tissue, y=percent_CD3, fill=treatment)) +
  geom_boxplot() + facet_wrap(~cell_type, scales = 'free') + theme(strip.text = element_text(size=6))

#####


CD3_tests <- CD3.long %>% 
  group_by(tissue, cell_type) %>%
  summarise(wilcox_p=wilcox.test(percent_CD3~treatment)$p.value, 
            t_p=t.test(percent_CD3~treatment)$p.value)

CD3_tests

###


###



# only cell types with a sig wilcox.test
CD3_sigs <- CD3_tests %>% filter(wilcox_p < 0.05)


CD3_sigs


###

CD3.long %>% filter(cell_type %in% CD3_sigs$cell_type) %>%
  ggplot(aes(x=tissue, y=percent_CD3, fill=treatment)) +
  geom_boxplot() + facet_wrap(~cell_type, scales = 'free') + theme(strip.text = element_text(size=6))

###


CD3_summary <- CD3.long %>% group_by(tissue, treatment, cell_type) %>% 
  summarise(avg=mean(percent_CD3), 
            std_dev= sd(percent_CD3), 
            num_obs= n(), 
            std_err=std_dev/sqrt(num_obs)) %>% filter(cell_type %in% CD3_sigs$cell_type)

CD3_summary


CD3_summary %>% ggplot(aes(x=tissue, y=avg, color=treatment)) +
  geom_point(size=3) + geom_errorbar(aes(ymin=avg-std_err, ymax=avg+std_err), width=.2) + 
  facet_wrap(~cell_type, scales = 'free') + theme(strip.text = element_text(size=6))



###


CD3_summary %>% filter(cell_type == 'CD4+CD8a-FoxP3+CD25+' & tissue == 'cec') %>% 
  ggplot(aes(x=treatment, y=avg, color=treatment)) + geom_errorbar(aes(ymin=avg-std_err, ymax= avg+std_err), width=0.2) +
  geom_point(size=5) + ggtitle('Percent CD3+CD4+CD8a-FoxP3+CD25+ cells') + ylab('Percent CD3+ cells')

###

CD3.wide <- CD3.long %>% mutate(set=paste(tissue, treatment,sep = '_')) %>% 
  spread(key = cell_type, value = percent_CD3)

CD3.wide


###


# you will need to run these commands if you want this next part to work: 

# install.packages('vegan') Awesome community ecology package
# devtools::install_github('jtrachsel/funfuns') 

library(funfuns)
NMDS <- NMDS_ellipse(metadata = CD3.wide[,1:4], OTU_table = CD3.wide[,-c(1:4)],grouping_set = 'set')

# how similar are these communities?
NMDS[[1]] %>% ggplot(aes(x=MDS1, y=MDS2, fill=set)) + geom_point(shape =21, size=2)



###


# join our sets with some segments...
NMDS[[1]] %>% ggplot(aes(x=MDS1, y=MDS2, fill=set)) + geom_point(shape =21, size=2) +
  geom_segment(aes(xend=centroidX, yend=centroidY, color=set))


###

# add in some standard error elipses...
NMDS[[1]] %>% ggplot(aes(x=MDS1, y=MDS2, fill=set)) + geom_point(shape =21, size=2) +
  geom_segment(aes(xend=centroidX, yend=centroidY, color=set)) +
  geom_path(data = NMDS[[2]], aes(x=NMDS1, y=NMDS2, color=group), inherit.aes = FALSE)

