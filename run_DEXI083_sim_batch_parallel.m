% Run DEXI_083 simulation pairs (9_10 to 15_16) in parallel.
% Each config is a single-interval run: first_time == last_time.

config_names = {
    'ours_DEXI083_sim_9_10'
    'ours_DEXI083_sim_10_11'
    'ours_DEXI083_sim_11_12'
    'ours_DEXI083_sim_12_13'
    'ours_DEXI083_sim_13_14'
    'ours_DEXI083_sim_14_15'
    'ours_DEXI083_sim_15_16'
};

nJobs = numel(config_names);

if nJobs <= 1
    fprintf('Single config detected; running serially without parallel pool.\n');
    run_one_config(config_names{1});
    fprintf('Single DEXI083 simulation config completed.\n');
    return;
end

pool = gcp('nocreate');
if isempty(pool)
    nWorkers = min(nJobs, feature('numcores'));
    fprintf('Starting parallel pool with %d workers...\n', nWorkers);
    pool = parpool(nWorkers);
else
    fprintf('Using existing pool with %d workers.\n', pool.NumWorkers);
end

status_ok = false(nJobs, 1);
err_msgs = strings(nJobs, 1);

parfor k = 1:nJobs
    cfg_name = config_names{k};
    try
        run_one_config(cfg_name);
        status_ok(k) = true;
    catch ME
        status_ok(k) = false;
        err_msgs(k) = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

for k = 1:nJobs
    if status_ok(k)
        fprintf('[OK] %s\n', config_names{k});
    else
        fprintf('[FAIL] %s\n%s\n', config_names{k}, err_msgs(k));
    end
end

fprintf('All DEXI083 simulation configs completed.\n');

function run_one_config(config_name)
    fprintf('Running config: %s\n', config_name);
    this_file = mfilename('fullpath');
    repo_root = fileparts(this_file);
    driver_path = fullfile(repo_root, 'driver_CAA.m');
    old_pwd = pwd;
    cleanup_obj = onCleanup(@() cd(old_pwd)); %#ok<NASGU>
    cd(repo_root);
    config_path = config_name; %#ok<NASGU>
    run(driver_path);
end
