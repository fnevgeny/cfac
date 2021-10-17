CREATE TEMPORARY VIEW _sessions_v AS
  SELECT sess.*, s.symbol, s.anum, s.mass,
    MIN(l.nele) AS nele_min, MAX(l.nele) AS nele_max
    FROM sessions AS sess
      INNER JOIN species AS s ON (sess.sid = s.sid)
      INNER JOIN levels AS l ON (s.sid = l.sid)
    GROUP BY sess.sid;

CREATE TEMPORARY VIEW _species_v AS
  SELECT s.sid, s.symbol, s.anum, s.mass,
         MIN(l.id) AS id_min, MAX(l.id) AS id_max
    FROM species AS s
      INNER JOIN levels AS l ON (s.sid = l.sid)
    GROUP BY s.sid;

CREATE TEMPORARY VIEW _cstates_v AS
  SELECT l.sid, l.nele, s.anum - l.nele + 1 AS zsp, MIN(e) AS e_gs,
         COUNT(id) AS nlevels
    FROM levels AS l
      INNER JOIN species AS s ON (l.sid = s.sid)
    GROUP BY s.sid, l.nele;

CREATE TEMPORARY VIEW _rtransitions_v AS
  SELECT rt.*,
         li.nele, lf.e - li.e AS de
    FROM rtransitions AS rt
      INNER JOIN levels AS li ON (rt.sid = li.sid AND rt.ini_id = li.id)
      INNER JOIN levels AS lf ON (rt.sid = lf.sid AND rt.fin_id = lf.id);

CREATE TEMPORARY VIEW _aitransitions_v AS
  SELECT at.sid, at.ini_id, at.fin_id, at.rate, li.nele
    FROM aitransitions AS at
      INNER JOIN levels AS li ON (at.sid = li.sid AND at.ini_id = li.id);

CREATE TEMPORARY VIEW _ctransitions_v AS
  SELECT t.sid, t.ini_id, t.fin_id, t.type, t.cid,
         li.nele AS ini_nele, lf.nele AS fin_nele, lf.e - li.e AS de
    FROM ctransitions AS t
      INNER JOIN levels AS li ON (t.sid = li.sid AND t.ini_id = li.id)
      INNER JOIN levels AS lf ON (t.sid = lf.sid AND t.fin_id = lf.id)
    GROUP BY t.cid;

CREATE TEMPORARY TABLE _cache_temp (
    cid      INTEGER NOT NULL REFERENCES ctransitions(cid) ON DELETE CASCADE,
    rate     REAL    NOT NULL
);

CREATE TEMPORARY VIEW _crates_cached_v AS
  SELECT ct.sid, ct.ini_id, ct.fin_id, ct.type,
         ini_nele, fin_nele, de, tt.rate
    FROM _ctransitions_v AS ct
    INNER JOIN _cache_temp AS tt ON (ct.cid = tt.cid);


CREATE TEMPORARY VIEW _states_v AS
  SELECT s.sid, s.fid, s.id, s.e, s.e - l.e AS de, s.mj,
         l.nele, l.name, l.g, l.vn, l.vl, l.p, l.ncomplex, l.sname, l.uta
    FROM states AS s
      INNER JOIN _levels_v AS l ON (s.sid = l.sid AND s.level_id = l.id);

CREATE TEMPORARY VIEW _rtransitions_m_v AS
  SELECT rt.*,
         si.nele, sf.e - si.e AS de
    FROM rtransitions_m AS rt
      INNER JOIN _states_v AS si ON (rt.sid = si.sid AND rt.fid = si.fid AND rt.ini_id = si.id)
      INNER JOIN _states_v AS sf ON (rt.sid = sf.sid AND rt.fid = sf.fid AND rt.fin_id = sf.id);
