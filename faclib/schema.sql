CREATE TABLE cfacdb (
    property TEXT UNIQUE NOT NULL,
    value INTEGER NOT NULL
);

CREATE TABLE fields (
    id      INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
    bf      REAL NOT NULL,
    ef      REAL NOT NULL,
    angle   REAL NOT NULL,
    CONSTRAINT unique_fields UNIQUE (bf, ef, angle)
);

INSERT INTO fields (id, bf, ef, angle) VALUES (0, 0.0, 0.0, 0.0);

CREATE TABLE sessions (
    sid     INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
    version INTEGER NOT NULL,
    tstamp  INTEGER NOT NULL,
    cmdline TEXT NOT NULL
);

CREATE TABLE species (
    sid    INTEGER NOT NULL REFERENCES sessions(sid) ON DELETE CASCADE,
    symbol TEXT    NOT NULL,
    anum   INTEGER NOT NULL,
    mass   REAL    NOT NULL
);

CREATE TABLE levels (
    sid      INTEGER NOT NULL REFERENCES sessions(sid) ON DELETE CASCADE,
    id       INTEGER NOT NULL,
    nele     INTEGER NOT NULL,
    name     TEXT    NOT NULL,
    e        REAL    NOT NULL,
    g        INTEGER NOT NULL,
    vn       INTEGER NOT NULL,
    vl       INTEGER NOT NULL,
    p        INTEGER NOT NULL,
    ibase    INTEGER,
    ncomplex TEXT    NOT NULL,
    sname    TEXT    NOT NULL,
    uta      BOOLEAN NOT NULL,
    PRIMARY KEY(sid, id),
    CONSTRAINT unique_name UNIQUE (sid, nele, name)
);

CREATE TABLE rtransitions (
    sid    INTEGER NOT NULL REFERENCES sessions(sid) ON DELETE CASCADE,
    ini_id INTEGER NOT NULL,
    fin_id INTEGER NOT NULL,
    mpole  INTEGER NOT NULL,
    rme    REAL    NOT NULL,
    mode   INTEGER NOT NULL,
    uta_de REAL,
    uta_sd REAL,
    FOREIGN KEY(sid, ini_id) REFERENCES levels(sid, id) ON DELETE CASCADE,
    FOREIGN KEY(sid, fin_id) REFERENCES levels(sid, id) ON DELETE CASCADE
);

CREATE TABLE aitransitions (
    sid    INTEGER NOT NULL REFERENCES sessions(sid) ON DELETE CASCADE,
    ini_id INTEGER NOT NULL,
    fin_id INTEGER NOT NULL,
    rate   REAL    NOT NULL,
    FOREIGN KEY(sid, ini_id) REFERENCES levels(sid, id) ON DELETE CASCADE,
    FOREIGN KEY(sid, fin_id) REFERENCES levels(sid, id) ON DELETE CASCADE
);

CREATE TABLE ctransitions (
    cid      INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
    sid      INTEGER NOT NULL REFERENCES sessions(sid) ON DELETE CASCADE,
    ini_id   INTEGER NOT NULL,
    fin_id   INTEGER NOT NULL,
    type     INTEGER NOT NULL,
    qk_mode  INTEGER NOT NULL,
    kl       INTEGER NOT NULL,
    ap0      REAL    NOT NULL,
    ap1      REAL    NOT NULL,
    ap2      REAL    NOT NULL,
    ap3      REAL    NOT NULL,
    FOREIGN KEY(sid, ini_id) REFERENCES levels(sid, id) ON DELETE CASCADE,
    FOREIGN KEY(sid, fin_id) REFERENCES levels(sid, id) ON DELETE CASCADE
);

CREATE TABLE cstrengths (
    cid      INTEGER NOT NULL REFERENCES ctransitions(cid) ON DELETE CASCADE,
    e        REAL    NOT NULL,
    strength REAL    NOT NULL
);


CREATE TABLE states (
    sid      INTEGER NOT NULL REFERENCES sessions(sid) ON DELETE CASCADE,
    fid      INTEGER NOT NULL REFERENCES fields(id) ON DELETE CASCADE,
    id       INTEGER NOT NULL,
    e        REAL    NOT NULL,
    level_id INTEGER NOT NULL,
    mj       INTEGER NOT NULL,
    PRIMARY KEY(sid, fid, id),
    FOREIGN KEY(sid, level_id) REFERENCES levels(sid, id) ON DELETE CASCADE
);

CREATE TABLE rtransitions_m (
    sid    INTEGER NOT NULL REFERENCES sessions(sid) ON DELETE CASCADE,
    fid    INTEGER NOT NULL REFERENCES fields(id) ON DELETE CASCADE,
    ini_id INTEGER NOT NULL,
    fin_id INTEGER NOT NULL,
    mpole  INTEGER NOT NULL,
    q      INTEGER NOT NULL,
    rme    REAL    NOT NULL,
    mode   INTEGER NOT NULL,
    FOREIGN KEY(sid, fid, ini_id) REFERENCES states(sid, fid, id) ON DELETE CASCADE,
    FOREIGN KEY(sid, fid, fin_id) REFERENCES states(sid, fid, id) ON DELETE CASCADE
);
