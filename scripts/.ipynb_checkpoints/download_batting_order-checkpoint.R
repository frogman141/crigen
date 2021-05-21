#!/home/baker02/miniconda3/envs/tidyscreen/bin/Rscript

# execute the analysis pipeline for the ER data
# each step will write intermediate output files to disk

#! How many cores per task?
#SBATCH --cpus-per-task=1
#! How much memory do you need?
#SBATCH --mem=4G
#SBATCH -J testing_dyngen_on_cluster
#! How much wallclock time will be required?
#SBATCH --time=1-00:00:00
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! General partition
#SBATCH -p general

library(progress)
library(textclean)
library(baseballr)
library(tidyverse)
library(argparse)

game_batting_order = function(game_id) {
  # get games batting order and remove first half of hyphened names
  batting_order = get_batting_orders(game_id) %>%
    rowwise %>%
    mutate(game_id = game_id,
           fullName = ifelse(length(unlist(strsplit(fullName, split = "-"))) == 2,
                             paste(unlist(strsplit(fullName, split = " "))[1],
                                   unlist(strsplit(unlist(strsplit(fullName, split = " "))[2], split="-"))[2]),
                             fullName))
  
  return (batting_order)
}

parser = ArgumentParser(description='Downloading Batting Order for all of the Games of a Given Year')
parser$add_argument('--start_date', help='Start Date for filtering games', type='character')
parser$add_argument('--end_date', help='End Date for filtering games', type='character')
arguments =  parser$parse_args()


end_date = as.Date(arguments$end_date)
start_date = as.Date(arguments$start_date)

game_ids = read.csv('/scratcha/fmlab/baker02/tidyscreen/data/baseball/2015_to_2020_game_ids.csv') 
game_ids = game_ids %>%
            mutate(game_date = as.Date(game_date)) %>% 
            filter(game_date >= start_date & game_date <= end_date)

message(paste0('Downloading Batting Orders from ', start_date, ' to ', end_date,'...'))

sleep_time = 1
batting_order = data.frame()
filename = paste0('games_batting_orders_', start_date,'_to_', end_date,'.csv')
data_dir = '/scratcha/fmlab/baker02/tidyscreen/data/baseball/batting_order'
filepath = file.path(data_dir, filename)

pb = progress_bar$new(total = length(game_ids$game_pk))

for (game_id in game_ids$game_pk) {
    bat_order = game_batting_order(game_id)
    batting_order = rbind(batting_order, bat_order)
    
    write.csv(batting_order, filepath, row.names=FALSE, append = T)
    
    Sys.sleep(sleep_time)
    pb$tick()
}

message(paste0('Saving Batting Orders from ', start_date, ' to ', end_date,'...'))

