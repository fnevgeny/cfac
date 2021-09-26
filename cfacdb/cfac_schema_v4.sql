CREATE TEMPORARY VIEW _levels_v AS
  SELECT l.sid, l.id, l.name, l.nele, cs.zsp, l.e, l.e - cs.e_gs AS e_rel,
         l.g, l.vn, l.vl, l.p, l.ncomplex, l.sname, l.uta
    FROM levels AS l
      INNER JOIN _cstates_v AS cs ON (l.sid = cs.sid AND l.nele = cs.nele);

