
#First step is to match the taxonomy used for the phylogeny to the taxonomy in the datafile. 
#First comparison indicated that there were 19 species which were in the database but not in the phylogeny, 
#which suggests there are mismatches in the taxonomy.

taxo<-read.csv("vertlife_taxonomies.csv")
data<-read.csv("mammal_data_imputed_withquality.csv")
taxo2<-subset(taxo, group=="Mammals")
taxo2$species<-sub(" ", "_", taxo2$scientificname)
data$species_name<-paste(data$genus, data$species, sep="_")
data<-data[, c("species_name", setdiff(names(data), "species_name"))]
setdiff(data$species_name, taxo2$species)

#Replace the species names from the database with those used in the phylogeny 
#(needs to be done just once as I will create a new dataframe with the correct names)   

#Dataframe with the old and new species names
lookup_df<-read.csv("old_new_species_names.csv")

# Match species names with the outdated names
match_idx<-match(data$species_name, lookup_df$species_name_anage)

#Replace species names where there is a match
data$species_name[!is.na(match_idx)]<-lookup_df$synonym_1[match_idx[!is.na(match_idx)]]

setdiff(data$species_name, taxo2$species)

write.csv(data, file="mammal_data_imputed_withqualityNames.csv")



#Now that the species names have been matched to the taxonomy of the phylogeny and adjusted in the databases 
#we can compare the phylogenies (downloaded from vertlife.org) to the databases and start the preliminary analyses. 

library(ape)
library(phangorn)
data.imp<-read.csv("mammal_data_imputed_withqualityNames.csv")
data.nas<-read.csv("mammal_data.csv")
phylo.all1<-read.nexus("MaxCredAll.nex")
phylo.dna1<-read.nexus("MaxCredDNA.nex")
setdiff(data.imp$species_name, phylo.all1$tip.label)
setdiff(phylo.all1$tip.label, data.imp$species_name)
#Drop the row for Balaena_mysticetus
data.imp<-data.imp[data.imp$species_name != "Balaena_mysticetus", ]
data.imp<-data.imp[data.imp$species != "wolfi", ]
```

#There are two species to be droped from the database, one is a whale (Balaena mysticetus) and the other 
#is a primate (Cercopithecus wolfi) because there is another species that matches the taxonomy in the tree. 

#All species in the dataset match the trees downloaded from the website. Now to check whether trees 
#estimated only for species for which DNA data is available match the species included in our database.

setdiff(data.imp$species_name, phylo.dna1$tip.label)

#There are 42 species for which we have data, at least in the file with imputed data, for which DNA 
#is not available for the phylogeny. Lets compare with the database for which data was not imputed droping all NAs.

lookup_df<-read.csv("old_new_species_names.csv")
data.nas$species_name<-paste(data.nas$genus, data.nas$species, sep="_")
data.nas<-data.nas[, c("species_name", setdiff(names(data.nas), "species_name"))]
data.nas2<-na.omit(data.nas)
#Also check whether the data frame without imputed data has incorrect names:
# Match species names with the outdated names
match_idx2<-match(data.nas2$species_name, lookup_df$species_name_anage)

#Replace species names where there is a match
data.nas2$species_name[!is.na(match_idx2)]<-lookup_df$synonym_1[match_idx2[!is.na(match_idx2)]]
setdiff(data.nas2$species_name, phylo.dna1$tip.label)
```

#When NAs are omitted from the database without imputation the phylogeny with only species for which 
#DNA data is available contains all species in that dataframe, after having adjusted the species names.
#All path analyses were run without imputed data and with the phylogeny with molecular data. 

#Can check whether the species included in the machine learning analyses and those in the phylogenetic path analysis 
#are distributed randomly or grouped in their respective phylogenies.   
#First the Machine Learning analyses:   

library(caper)

# Assuming species_long is the longer vector and species_short is the shorter one
species_long <- phylo.dna1$tip.label
species_short <- data.nas2$species_name

# Create the match indicator: 1 if species in species_long is in species_short, 0 otherwise
match_indicator <- as.integer(species_long %in% species_short)

# Create the data frame
result_df <- data.frame(Species = species_long, Match = match_indicator)

# View the result
head(result_df)

#Remove node labels from the DNA tree:
phylo.dna1$node.label <- NULL

phylo.d(data=result_df, phy=phylo.dna1, names.col=Species, binvar=Match, permut = 1000, rnd.bias=NULL)


#### 1. Test of the effects of data quality and sample size on estimates of lifespan   

#Test whether data quality and sample size (provided in AnAge as categories) has any 
#association with lifespan estimates of species. I will also test the association between body size 
#and sample size, just to check whether large species, possibly more studied, have larger sample sizes.

t1<-aov(log10(adult_mass_g)~sample.size, data=data.nas2)
summary(t1)
TukeyHSD(t1)
t2<-aov(log10(max_longevity_d)~sample.size, data = data.nas2)
summary(t2)
TukeyHSD(t2)
t3<-aov(log10(max_longevity_d)~sample.size+log10(adult_mass_g), data = data.nas2)
summary(t3)
TukeyHSD(t3, "sample.size")
#These models however do not inform me regarding the magnitude of the effect of sample size on estimates of maximum longevity:
library(phylolm)
library(ggplot2)
to.drop<-setdiff(phylo.dna1$tip.label, data.nas2$species_name)
phylo.dna.cut<-drop.tip(phylo.dna1, to.drop)
setdiff(phylo.dna.cut$tip.label, data.nas2$species_name)
row.names(data.nas2)<-data.nas2$species_name
t4<-phylolm(log10(adult_mass_g)~sample.size, data = data.nas2, phy = phylo.dna.cut, model="lambda")
summary(t4)
t5<-phylolm(log10(max_longevity_d)~sample.size, data = data.nas2, phy = phylo.dna.cut, model="lambda")
summary(t5)
t6<-phylolm(log10(max_longevity_d)~sample.size+log10(adult_mass_g), data = data.nas2, phy = phylo.dna.cut, model="lambda")
summary(t6)
t7<-phylolm(log10(max_longevity_d)~log10(adult_mass_g), data = data.nas2, phy = phylo.dna.cut, model="lambda")
summary(t7)
#The difference in R2 between models including and not including sample size is very small = 0.0518

plot(t7$residuals~as.factor(data.nas2$sample.size))

#To export to PDF
#pdf("violin_plot_residuals.pdf", width = 6, height = 5)

ggplot(data.nas2, aes(x = as.factor(data.nas2$sample.size), y = t7$residuals)) +
  geom_violin(fill = "lightgreen", trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  labs(x = "Sample Size Category", y = "Model Residuals",
       title = "Distribution of Residuals by Sample Size Category") +
  theme_minimal()

#dev.off()

df <- data.frame(
  residuals = t7$residuals,
  sample_size_category = as.factor(data.nas2$sample.size)
)
agg_results <- aggregate(residuals ~ sample_size_category, data = df,
                         FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_df <- data.frame(
  sample_size_category = agg_results$sample_size_category,
  mean_residual = agg_results$residuals[, "mean"],
  sd_residual = agg_results$residuals[, "sd"]
)

print(agg_df)

anova_model <- aov(t7$residuals ~ as.factor(data.nas2$sample.size), data = data.nas2)
summary(anova_model)

#Test whether the variance of model residuals differs across sample size categories — which would indicate heteroskedasticity
library(car)
leveneTest(t7$residuals ~ as.factor(data.nas2$sample.size), data = data.nas2)
```

###Phylogenetic path analysis

#First define a set of models, the first group will involve models where alometry (i.e. body mass) 
#determines the life-history traits, and then variations of these life-history traits, as well as 
#dispersal and body mass itself, will affect longevity. The second group of models will involve life-history traits 
#determined by brain size (reflecting the idea of energetic constraints imposed by brain size 
#(see e.g. Barton and Capellini 2011 and Gonzalez-Voyer et al 2016) and then these life-history traits, 
#and brain size as well as dispersal, affecting longevity. Finally, the last set of models will involve a combination 
#of allometry and brain size influencing life-history traits. Such a combination of models will allow us to both distinguish 
#between purely allometric, purely energetic cost of brain size effects on life-history traits, as well as to analyze the causal 
#effects of body size, brain size, life-history and dispersal on longevity. 

#Create another database with simplified variable names, as I will need to type these in multiple times for the path models.

#Take a two-step approach for the path models, first to define the relationships between body mass, brain size and the life-history traits 
#and then fix those to the ones found to fit the data best - and for which all conditional independencies are met - and then 
#explore the relationships with longevity. So below is a set of path models describing potential causal relationships between mass, brain size and life-history traits. 

row.names(data.nas2)<-data.nas2$species_name
setdiff(data.nas2$species_name, phylo.dna1$tip.label)
data.nas3<-subset(data.nas2, select=c("species_name", "adult_mass_g", "adult_brain_mass_g","gestation_length_d","litters_per_year_n", "weaning_age_d", "female_maturity_d", "max_longevity_d"))
colnames(data.nas3)<-c("species_name", "mass", "brain", "gestation", "litters", "weaning", "fmaturity", "longevity")
data.nas3[2:8]<-lapply(data.nas3[2:8], log10)
row.names(data.nas3)<-data.nas3$species_name

#Need to crop the phylogeny to match the species for which we have data:   


to.dropDNA<-setdiff(phylo.dna1$tip.label, row.names(data.nas3))
phylo.dna2<-drop.tip(phylo.dna1, to.dropDNA)
setdiff(phylo.dna2$tip.label, row.names(data.nas3))

library(phylopath)
models.lh<-define_model_set(
  one = c(weaning~mass, fmaturity~weaning+gestation, gestation~mass, litters~gestation),
  two = c(weaning~mass+fmaturity, fmaturity~mass, gestation~mass, litters~mass),
  three = c(litters~mass+fmaturity, weaning~mass+fmaturity,fmaturity~mass, gestation~fmaturity+weaning),
  four = c(weaning~brain, fmaturity~weaning+gestation, gestation~brain, litters~gestation),
  five = c(weaning~brain+fmaturity, fmaturity~brain, gestation~brain, litters~brain),
  six = c(litters~mass+fmaturity, weaning~mass+fmaturity,fmaturity~mass, gestation~fmaturity+brain),
  seven = c(litters~mass+weaning+fmaturity, weaning~brain+fmaturity,fmaturity~mass, gestation~fmaturity+brain),
  eight = c(litters~mass+fmaturity, weaning~brain+fmaturity+gestation,fmaturity~mass, gestation~fmaturity+brain),
  nine = c(litters~mass+weaning+fmaturity, weaning~brain+fmaturity+gestation,fmaturity~mass, gestation~fmaturity+brain),
  ten = c(litters~mass+weaning+gestation, weaning~brain+fmaturity+gestation,fmaturity~brain, gestation~fmaturity+brain),
  eleven = c(litters~mass, weaning~brain+mass+gestation, fmaturity~mass, gestation~fmaturity+brain),
  twelve = c(litters~mass, weaning~brain+mass+gestation, fmaturity~mass, gestation~brain),
  thirteen = c(litters~mass, weaning~brain+mass+gestation, fmaturity~weaning+mass, gestation~brain),
  fourteen = c(litters~mass, weaning~brain+mass+gestation, fmaturity~weaning+gestation, gestation~brain),
  fifteen = c(litters~mass+weaning+fmaturity, weaning~brain+gestation+mass,fmaturity~mass+gestation, gestation~brain),
  sixteen = c(litters~mass+weaning+fmaturity, weaning~brain+gestation+mass, fmaturity~mass+brain+weaning, gestation~brain),
  seventeen = c(litters~mass+weaning+fmaturity, weaning~brain+gestation+mass, fmaturity~mass+brain+weaning+gestation, gestation~brain),
  eighteen = c(litters~mass+weaning+fmaturity, weaning~brain+mass+gestation,fmaturity~mass+brain, gestation~fmaturity+brain),
  nineteen = c(litters~mass+weaning+fmaturity, weaning~brain+gestation+mass, fmaturity~mass+brain+weaning+gestation, gestation~brain+mass),
  twenty = c(litters~mass+weaning+fmaturity, weaning~brain+gestation+mass, fmaturity~brain+weaning+gestation, gestation~brain+mass),
    .common = c(brain~mass)
)

#Determine what is the best model for the relationships among life-history traits. 

result.lh.DNA <- phylo_path(models.lh, data = data.nas3, tree = phylo.dna2, model = 'lambda')


#To check conditional independencies (e.g.):   

result.lh.DNA$d_sep$seventeen
result.lh.DNA$d_sep$nineteen
result.lh.DNA$d_sep$twenty
result.lh.DNA$d_sep$eighteen
```


(s.DNA<-summary(result.lh.DNA))

#Four models can be considered potential causal models. 

(best.lh.DNA<-best(result.lh.DNA))

#Can look at the different models that fulfilled the conditional independencies to see how they differ. 
#Need to fit each of them again and then we can plot them.

d_fit_lhDNA17<-est_DAG(models.lh$seventeen, data = data.nas3, tree = phylo.dna2, model = 'lambda')

d_fit_lhDNA19<-est_DAG(models.lh$nineteen, data = data.nas3, tree = phylo.dna2, model = 'lambda')

d_fit_lhDNA20<-est_DAG(models.lh$twenty, data = data.nas3, tree = phylo.dna2, model = 'lambda')

d_fit_lhDNA18<-est_DAG(models.lh$eighteen, data = data.nas3, tree = phylo.dna2, model = 'lambda')

#Can now plot each model:   

plot(d_fit_lhDNA17, text_size=3, box_x=13, box_y=10, curvature=0.15)

plot(d_fit_lhDNA19, text_size=3, box_x=13, box_y=10, curvature=0.15)

plot(d_fit_lhDNA20, text_size=3, box_x=13, box_y=10, curvature=0.15)

plot(d_fit_lhDNA18, text_size=3, box_x=13, box_y=10, curvature=0.15)

#Four models which fulfill the conditional independencies, and therefore can be considered possible 
#causal models. The best model is #17, the second-best model is #19, with a deltaCIC of 1.92, followed 
#by models 20 and 18 with deltaCIC of > 2. We can use model 17 for the path analyses of 
#the relationship between body size, brain size and life history traits with lifespan. 


```{r echo=TRUE}
models.fullDNA<-define_model_set(
  one = c(longevity~fmaturity),
  two = c(longevity~fmaturity+mass),
  three = c(longevity~fmaturity+brain),
  four = c(longevity~weaning),
  five = c(longevity~weaning+mass),
  six = c(longevity~weaning+brain),
  seven = c(longevity~gestation),
  eight = c(longevity~gestation+mass),
  nine = c(longevity~gestation+brain),
  ten = c(longevity~litters),
  eleven = c(longevity~litters+mass),
  twelve = c(longevity~litters+brain),
  thirteen = c(longevity~fmaturity+weaning),
  fourteen = c(longevity~fmaturity+weaning+mass),
  fifteen = c(longevity~fmaturity+weaning+brain),
  sixteen = c(longevity~weaning+gestation),
  seventeen = c(longevity~weaning+gestation+mass),
  eighteen = c(longevity~weaning+gestation+brain),
  nineteen = c(longevity~gestation+litters),
  twenty = c(longevity~gestation+litters+mass),
  twentyone = c(longevity~gestation+litters+brain),
  twentytwo = c(longevity~fmaturity+gestation),
  twentythree = c(longevity~fmaturity+gestation+mass),
  twentyfour = c(longevity~fmaturity+gestation+brain),
  .common = c(brain~mass, litters~mass+weaning+fmaturity, weaning~brain+gestation+mass, fmaturity~mass+brain+weaning+gestation, gestation~brain)
)

#Run the models

result.full.DNA <- phylo_path(models.fullDNA, data = data.nas3, tree = phylo.dna2, model = 'lambda')


#Get summary of the results:   

s.DNA<-summary(result.full.DNA)
s.DNA

#Check the best model:   

(best_model<-best(result.full.DNA))

#Organize the variables to make plotting the models easier

m1<-data.frame(
  name=c("mass", "brain", "weaning", "gestation", "litters", "fmaturity", "dispersal", "longevity"),
  x=c(2,4,3,4,2,5,1,3),
  y=c(3,3,2,2,2,2,2,1)
)

#Plot the best model:

plot(best_model, text_size=3, box_x=13, box_y=10, curvature=0.15, manual_layout = m1)

#Check the averaged model with a cut-off of delta CIC of 2. 

average.mod.DNA<-average(result.full.DNA)
average.mod.DNA

plot(average.mod.DNA, text_size=3, box_x=13, box_y=10, curvature=0.15, manual_layout = m1)

#It is also worth looking at each of the specific models that compose the averaged model (these are model 21 and model 24):   

d_fit21<-est_DAG(models.fullDNA$twentyone, data = data.nas3, tree = phylo.dna2, "lambda", boot=100)
d_fit24<-est_DAG(models.fullDNA$twentyfour, data = data.nas3, tree = phylo.dna2, "lambda", boot=100)

#Can plot model 24 (since model 21 is the best model and we already plotted it):   

plot(d_fit24, text_size=3, box_x=13, box_y=10, curvature=0.15, manual_layout = m1)

#Export the table of results to an editable Word document:   

library(flextable)
library(officer)

# Get model comparison table
model_table <- summary(result.full.DNA)

is_numeric <- sapply(model_table, is.numeric)
model_table[is_numeric] <- lapply(model_table[is_numeric], function(x) round(x, 3))

# Clean up column names for clarity (optional)
names(model_table) <- gsub("\\.", " ", names(model_table))

# Create a flextable object
ft <- flextable(model_table)
ft <- autofit(ft)

# Export to Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "phylopath_model_comparison_fmat.docx")

#Repeat for the table of the life history models:   

# Get model comparison table
model_table_lh <- summary(result.lh.DNA)

is_numeric <- sapply(model_table_lh, is.numeric)
model_table_lh[is_numeric] <- lapply(model_table_lh[is_numeric], function(x) round(x, 3))

# Clean up column names for clarity (optional)
names(model_table_lh) <- gsub("\\.", " ", names(model_table_lh))

# Create a flextable object
ft_lh <- flextable(model_table_lh)
ft_lh <- autofit(ft_lh)

# Export to Word
doc_lh <- read_docx()
doc_lh <- body_add_flextable(doc_lh, ft_lh)
print(doc_lh, target = "phylopath_model_comparison_lh_fmat.docx")
