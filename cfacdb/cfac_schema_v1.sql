CREATE TEMPORARY VIEW _cstrengths_v AS
  SELECT ct.sid, ct.ini_id, ct.fin_id, ct.type,
         cs.e, cs.strength,
         li.nele AS ini_nele, lf.nele AS fin_nele,
         lf.e - li.e AS de, ct.ap0, ct.ap1
    FROM cstrengths AS cs
      INNER JOIN ctransitions AS ct ON (cs.cid = ct.cid)
      INNER JOIN levels AS li ON (ct.sid = li.sid AND ct.ini_id = li.id)
      INNER JOIN levels AS lf ON (ct.sid = lf.sid AND ct.fin_id = lf.id)
    WHERE cs.strength > 0;
