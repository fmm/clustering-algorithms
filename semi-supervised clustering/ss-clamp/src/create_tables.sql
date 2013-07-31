CREATE TABLE IF NOT EXISTS algorithm (
    -- extracted from config file
    sha1 TEXT,
    info TEXT,
    seed INTEGER,
    time_limit REAL,
    class_variable INTEGER,
    alpha REAL,
    label_percentage REAL,
    initializations INTEGER,
    clusters INTEGER,
    individuals INTEGER,
    prototypes INTEGER,
    eps REAL,
    relevance_v REAL,
    iterations INTEGER,
    self_training_in REAL,
    self_training_out REAL,
    -- generated at the end
    used_pwc_file TEXT,
    used_clusters INTEGER,
    used_prototypes INTEGER,
    used_alpha REAL,
    best_initialization INTEGER,
    PRIMARY KEY(sha1)
    );

CREATE TABLE IF NOT EXISTS input (
    algorithm_id TEXT,
    file TEXT,
    PRIMARY KEY(algorithm_id, file),
    FOREIGN KEY(algorithm_id) REFERENCES algorithm(sha1)
    );

CREATE TABLE IF NOT EXISTS answer (
    algorithm_id TEXT,
    initialization INTEGER,
    iteration INTEGER,
    criterion REAL,
    restriction REAL,
    accuracy REAL,
    adjusted_rand_index REAL,
    f_measure REAL,
    qr REAL,
    PRIMARY KEY(algorithm_id,initialization,iteration),
    FOREIGN KEY(algorithm_id) REFERENCES algorithm(sha1)
    );

CREATE TABLE IF NOT EXISTS partition (
    algorithm_id TEXT,
    individual INTEGER,
    cluster INTEGER,
    value REAL,
    PRIMARY KEY(algorithm_id,individual,cluster),
    FOREIGN KEY(algorithm_id) REFERENCES algorithm(sha1)
    );

CREATE TABLE IF NOT EXISTS relevance (
    algorithm_id TEXT,
    cluster INTEGER,
    matrix INTEGER,
    value REAL,
    PRIMARY KEY(algorithm_id,cluster,matrix),
    FOREIGN KEY(algorithm_id) REFERENCES algorithm(sha1)
    );
