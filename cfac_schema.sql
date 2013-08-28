CREATE TEMPORARY VIEW _species_v AS
  SELECT s.sid, s.symbol, s.anum, s.mass,
         MIN(l.id) AS id_min, MAX(l.id) AS id_max
    FROM species AS s
      INNER JOIN levels AS l ON (s.sid = l.sid);

CREATE TEMPORARY VIEW _cstates_v AS
  SELECT l.sid, l.nele, s.anum - l.nele + 1 AS zsp, MIN(e) AS e_gs,
         COUNT(id) AS nlevels
    FROM levels AS l
      INNER JOIN species AS s ON (l.sid = s.sid)
    GROUP BY l.nele;

CREATE TEMPORARY VIEW _levels_v AS
  SELECT l.sid, l.id, l.name, l.nele, cs.zsp, l.e - cs.e_gs AS e,
         l.g, l.vn, l.vl, l.p, l.ncomplex, l.sname
    FROM levels AS l
      INNER JOIN _cstates_v AS cs ON (l.sid = cs.sid AND l.nele = cs.nele);

CREATE TEMPORARY VIEW _rtransitions_v AS
  SELECT rt.sid, rt.ini_id, rt.fin_id, rt.mpole, rt.rme,
         li.nele, li.e - lf.e AS de
    FROM rtransitions AS rt
      INNER JOIN levels AS li ON (rt.sid = li.sid AND rt.ini_id = li.id)
      INNER JOIN levels AS lf ON (rt.sid = lf.sid AND rt.fin_id = lf.id);

CREATE TEMPORARY VIEW _aitransitions_v AS
  SELECT at.sid, at.ini_id, at.fin_id, at.rate, li.nele
    FROM aitransitions AS at
      INNER JOIN levels AS li ON (at.sid = li.sid AND at.ini_id = li.id);

CREATE TEMPORARY VIEW _ctransitions_v AS
  SELECT t.sid, t.ini_id, t.fin_id, t.type, t.cid, COUNT(s.cid) AS ndata,
         t.ap0, t.ap1,
         li.nele AS ini_nele, lf.nele AS fin_nele, lf.e - li.e AS de
    FROM ctransitions AS t
      INNER JOIN cstrengths AS s ON (s.cid = t.cid)
      INNER JOIN levels AS li ON (t.sid = li.sid AND t.ini_id = li.id)
      INNER JOIN levels AS lf ON (t.sid = lf.sid AND t.fin_id = lf.id)
    WHERE s.strength > 0
    GROUP BY t.cid;

CREATE TEMPORARY VIEW _cstrengths_v AS
  SELECT ct.sid, ct.ini_id, ct.fin_id, ct.type,
         cs.e, cs.strength,
         li.nele AS ini_nele, lf.nele AS fin_nele,
         li.e - lf.e AS de, ct.ap0, ct.ap1
    FROM cstrengths AS cs
      INNER JOIN ctransitions AS ct ON (cs.cid = ct.cid)
      INNER JOIN levels AS li ON (ct.sid = li.sid AND ct.ini_id = li.id)
      INNER JOIN levels AS lf ON (ct.sid = lf.sid AND ct.fin_id = lf.id)
    WHERE cs.strength > 0;
