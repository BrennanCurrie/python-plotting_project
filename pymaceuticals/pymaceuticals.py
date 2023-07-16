# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from scipy.stats import linregress

# Study data files
mouse_metadata_path = "/Users/brennancurrie/Desktop/My_Code/Class_Materials/Challenges/Module_5_Challenge/pymaceuticals/data/mouse_metadata.csv"
study_results_path = "/Users/brennancurrie/Desktop/My_Code/Class_Materials/Challenges/Module_5_Challenge/pymaceuticals/data/study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
datacombined = pd.merge(study_results, mouse_metadata, how="left", on=["Mouse ID"])

# Display the data table for preview
datacombined

# Checking the number of mice.
micecount = len(datacombined['Mouse ID'].unique())
micecount

# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
grp = datacombined.loc[datacombined.duplicated(subset=['Mouse ID', 'Timepoint']),'Mouse ID']
duplicate_mouse_ids = grp.unique()
duplicate_mouse_ids

# Optional: Get all the data for the duplicate mouse ID. 
dupegrp = datacombined.loc[(datacombined["Mouse ID"] == "g989")]

dupegrp

# Create a clean DataFrame by dropping the duplicate mouse by its ID.
cleandataframe = datacombined.loc[(datacombined["Mouse ID"] != "g989")]
cleandataframe

# Checking the number of mice in the clean DataFrame.
micecountclean = len(cleandataframe['Mouse ID'].unique())
micecountclean

#SUMMARY STATISTICS
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.

tvmean = cleandataframe.groupby('Drug Regimen').mean()['Tumor Volume (mm3)']
tvmedian = cleandataframe.groupby('Drug Regimen').median()['Tumor Volume (mm3)']
tvvariance = cleandataframe.groupby('Drug Regimen').var()['Tumor Volume (mm3)']
tvsd = cleandataframe.groupby('Drug Regimen').std()['Tumor Volume (mm3)']
tvsem = cleandataframe.groupby('Drug Regimen').sem()['Tumor Volume (mm3)']

sumstatdf = pd.DataFrame({"Mean Tumor Volume":tvmean,
                          "Median Tumor Volume": tvmedian,
                          "Tumor Volume Variance": tvvariance,
                          "Tumor Volume Std. Dev.": tvsd,
                          "Tumor Volume Std. Error": tvsem
                         })
sumstatdf

#BONUS/OPTIONAL
# Generate a summary statistics table of mean, median, variance, standard deviation, 
# and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line.

druggrp = cleandataframe.groupby('Drug Regimen')['Tumor Volume (mm3)']
druggrp.agg(['mean','median', 'var','std', 'sem'])

# BAR & PIE CHARTS
# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.

miceperdrug = cleandataframe['Drug Regimen'].value_counts()

miceperdrug.plot(kind="bar",figsize=(5,2))

plt.xlabel('Drug Regimen')
plt.ylabel('Number Of Mice Tested')

# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.

xvalues = miceperdrug.index.values
yvalues = miceperdrug.values

plt.bar(xvalues, yvalues)
plt.xlabel('Drug Regimen')
plt.ylabel('Number Of Mice Tested')
plt.xticks(rotation=90)
plt.show()

# Generate a pie plot showing the distribution of female versus male mice using pyplot
gendercounts = cleandataframe['Sex'].value_counts()
labels = ['Male', 'Female']
explode = (0,0)
plt.pie(gendercounts,explode=explode,labels=labels,autopct='%1.1f%%')
plt.ylabel('Sex')

# Generate a pie plot showing the distribution of female versus male mice using Pandas

gendercounts.plot(kind = 'pie', autopct='%1.1f%%')


#QUARTILES, OUTLIERS, AND BOXPLOTS

# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
#Start by getting the last (greatest) timepoint for each mouse

lasttime = cleandataframe.groupby('Mouse ID').max()['Timepoint']
lasttimedf = pd.DataFrame(lasttime)
lastdf = lasttimedf.reset_index()                     
lastdf

#Merge this group df with the original DataFrame to get the tumor volume at the last timepoint

newdf = pd.merge(cleandataframe, lastdf, how="right", on=["Mouse ID", 'Timepoint'])
newdf

# Put treatments into a list for for loop (and later for plot labels)
druglist = ['Capomulin','Ramicane','Infubinol','Ceftamin']

# Create empty list to fill with tumor vol data (for plotting)
emptylist = []

#Locate the rows which contain mice on each drug and get the tumor volumes 
#looping through your "druglist" varaible, selecting the drug regimen and its corresponding volume from "newdf"
#using iloc, and appending it to the empty list. 

for drug in druglist:
    final_tumor_vol = newdf.loc[newdf["Drug Regimen"] == drug, 'Tumor Volume (mm3)']
    
    # add subset 
    emptylist.append(final_tumor_vol)

    # Calculate the IQR and quantitatively determine if there are any potential outliers. 
    tumorvols = newdf['Tumor Volume (mm3)']
    tumorvolsqt = tumorvols.quantile([.25, .5, .75])

    upperq = tumorvolsqt[0.75] 
    lowerq = tumorvolsqt[0.25]
    innerqr = upperq-lowerq 

    outerbounds = upperq + (1.5 * innerqr)
    innerbounds = lowerq - (1.5 * innerqr)

    # Determine outliers using upper and lower bounds

    outliers = final_tumor_vol.loc[(final_tumor_vol < innerbounds) | (final_tumor_vol > outerbounds)]
    print(f"{drug}'s potential outliers: {outliers}")


# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.

fig1, ax1 = plt.subplots()
ax1.boxplot(emptylist)
ax1.set_ylabel('Tumor Volume (mm3)')
plt.xticks([1, 2, 3, 4], druglist)
plt.show()

#LINE AND SCATTER PLOTS
#Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin l509
#Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen

capomulin = cleandataframe.loc[(cleandataframe["Drug Regimen"] == "Capomulin")]
mouseid = capomulin.loc[(capomulin["Mouse ID"] == "l509")]

liney = mouseid['Tumor Volume (mm3)']
linex = mouseid['Timepoint']

plt.plot(linex,liney)

plt.xlabel('Timepoint')
plt.ylabel('Tumor Volume (mm3)')
plt.title('Capomulin treatment of mouse l509')



#HINT
#use the capomulin dataframe and use group by and mean()to calculate the averages, you will group by mouse id.
#store this in a new dataframe and name is accordingly and then you can use the .scatter function 
#to plot the weight against the tumor volume

scaty = capomulin.groupby('Mouse ID').mean()['Tumor Volume (mm3)']
scatx = capomulin.groupby('Mouse ID').mean()['Weight (g)']

plt.scatter(x=scatx, y=scaty)
plt.ylabel('Tumor Volume (mm3)')
plt.xlabel('Weight (g)')

#CORRELATION & REGRESSION
# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen

correlation = st.pearsonr(scatx,scaty)

(slope, intercept, rvalue, pvalue, stderr) = linregress(scatx,scaty)
regressvalues = scatx * slope + intercept 

print(f'The correlation between mouse weight and the average tumor volume is {round(correlation[0],2)}')

#plot scatter
plt.scatter(x=scatx, y=scaty)
plt.ylabel('Tumor Volume (mm3)')
plt.xlabel('Weight (g)')


#plot corr/reg line
plt.plot(scatx,regressvalues,'r-')

plt.show()