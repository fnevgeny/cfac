CREATE TABLE IF NOT EXISTS crates (
    cid      INTEGER NOT NULL,
    t        REAL    NOT NULL,
    rate     REAL    NOT NULL,
    UNIQUE (cid, t)
);
