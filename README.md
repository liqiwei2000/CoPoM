# CoPoM
The program in written in R. You need R to run the code. Library mcclust and lattice are required to get and show the PPM estimate of z.

1. how to load the data

One simulated data set and the real data used in the paper are ready to use. They are simulated_data.RData and real_data.RData, respectively. To generate different simulated data, please change the settings in the code and run simulated_data_generator.R. To load your own high-dimensional count data, please prepare the data in a n-by-p count matrix, where n is the number of rows (observations) and p is the number of columns (features).

2. how to run the code

To run the code, please compile main.R. First, the users need to define the hyper parameters, environments, and algorithm settings (such as iterations) by themselves, although we provide the default values. Please note the notations are as the same as those in the papers.

3. how to analyze the result

To analyze the result, please compile evaluation.R. You can get both numerical and visualized results for the parameters of interests.
