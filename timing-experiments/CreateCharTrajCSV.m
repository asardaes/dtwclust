# Download dataset
urlwrite("https://archive.ics.uci.edu/ml/machine-learning-databases/character-trajectories/mixoutALL_shifted.mat",
         "CharTraj.mat")
# Load it, each element of the "mixout" cell array is a multivariate time series
load CharTraj.mat;
# Make directory if it doesn't exist
dir_name = "CharTrajCSV";
if (~isdir(dir_name))
  mkdir(dir_name);
endif

# Function to write one series to the CSV
function append_series(fid, series)
  # For whatever reason, I originally used 5 significant digits when saving.
  # Remove the last two parameters if you want full precision.
  dlmwrite(fid, series, "delimiter", ",", "append", "on", "precision", 5);
endfunction

# Function to map the variables in one multivariate series to the above
function write_multivariate_series(series, fid1, fid2, fid3)
  append_series(fid1, series(1,:));
  append_series(fid2, series(2,:));
  append_series(fid3, series(3,:));
endfunction

# Loop across labels
for (id = unique(consts.charlabels))
  # Which series correspond to current id
  indices = consts.charlabels == id;
  # Which character is it
  character = toupper(consts.key{id});
  # Subset of cell array
  series = mixout(indices);
  # Will write each variable in a different CSV file
  fid_velX = fopen([dir_name, "/", character, "tipVelX.csv"], "w");
  fid_velY = fopen([dir_name, "/", character, "tipVelY.csv"], "w");
  fid_force = fopen([dir_name, "/", character, "tipForce.csv"], "w");
  # Write series
  cellfun(@write_multivariate_series, 
          series, 
          {fid_velX},
          {fid_velY},
          {fid_force});
  # Close files
  fclose(fid_velX);
  fclose(fid_velY);
  fclose(fid_force);
endfor
