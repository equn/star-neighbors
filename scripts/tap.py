"""Minimal TAP (Table Access Protocol) client for the Gaia and SIMBAD archives.

Uses only the stdlib so the pipeline has no install step. Sync queries for small
result sets, async jobs for anything large (Gaia caps sync at ~2000 rows).
"""

from __future__ import annotations

import csv
import http.client
import io
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap"
SIMBAD_TAP = "https://simbad.u-strasbg.fr/simbad/sim-tap"
# VizieR serves the published encounter catalogs this project validates
# against. Table names there contain '/' and '+' and must be double-quoted in
# ADQL: SELECT * FROM "J/A+A/616/A37/table23".
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"

UA = {"User-Agent": "star-neighbors/1.0 (research; encounter analysis)"}


# The ESA and SIMBAD endpoints intermittently reset connections under load;
# a long fetch that dies two thirds of the way through is expensive, so retry.
RETRIES = 6
BACKOFF = 5.0

# The failure that actually bit: ESA closed a chunked response half-delivered,
# which surfaces as http.client.IncompleteRead - an HTTPException, and so not
# caught by the URLError/ConnectionError/TimeoutError triple this used to list.
# Two retries were consumed by read timeouts and the third died on the truncated
# body, losing a completed server-side job. Catch the whole family.
TRANSIENT = (urllib.error.URLError, ConnectionError, TimeoutError,
             http.client.HTTPException, OSError)


def _with_retry(fn, what: str):
    for attempt in range(RETRIES):
        try:
            return fn()
        except TRANSIENT as exc:
            if attempt == RETRIES - 1:
                raise
            wait = BACKOFF * (attempt + 1)
            print(f"    {what} failed ({type(exc).__name__}: {exc}); "
                  f"retrying in {wait:.0f}s")
            time.sleep(wait)


def _post(url: str, data: dict, timeout: int = 900) -> bytes:
    def go():
        req = urllib.request.Request(
            url, data=urllib.parse.urlencode(data).encode(), headers=UA
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read()
    return _with_retry(go, "POST")


def _get(url: str, timeout: int = 900) -> bytes:
    def go():
        with urllib.request.urlopen(
            urllib.request.Request(url, headers=UA), timeout=timeout
        ) as resp:
            return resp.read()
    return _with_retry(go, "GET")


def sync_query(query: str, base: str = GAIA_TAP, timeout: int = 900) -> list[dict]:
    """Run a synchronous ADQL query, returning rows as dicts of strings."""
    raw = _post(
        base + "/sync",
        {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "csv",
            "QUERY": query,
        },
        timeout=timeout,
    )
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8", "replace"))))


def async_query(
    query: str, base: str = GAIA_TAP, poll: float = 3.0, max_wait: float = 3600
) -> list[dict]:
    """Submit an async ADQL job, poll to completion, return rows as dicts."""
    raw = _post(
        base + "/async",
        {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "csv",
            "PHASE": "RUN",
            "QUERY": query,
        },
    )
    root = ET.fromstring(raw)
    ns = {"uws": "http://www.ivoa.net/xml/UWS/v1.0"}
    job_id = root.find("uws:jobId", ns).text
    job_url = f"{base}/async/{job_id}"

    started = time.time()
    while True:
        phase = _get(job_url + "/phase").decode().strip()
        if phase in ("COMPLETED", "ERROR", "ABORTED"):
            break
        if time.time() - started > max_wait:
            raise TimeoutError(f"TAP job {job_id} still {phase} after {max_wait}s")
        time.sleep(poll)

    if phase != "COMPLETED":
        raise RuntimeError(f"TAP job {job_id} ended in {phase}:\n{_get(job_url).decode()[:4000]}")

    body = _get(job_url + "/results/result").decode("utf-8", "replace")
    return list(csv.DictReader(io.StringIO(body)))


def f(row: dict, key: str) -> float | None:
    """Parse a possibly-empty CSV cell as float."""
    v = row.get(key)
    if v is None or v == "" or v.lower() in ("nan", "null", "none"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


if __name__ == "__main__":
    rows = sync_query("SELECT TOP 3 source_id, parallax, radial_velocity "
                      "FROM gaiadr3.gaia_source WHERE parallax > 500")
    for r in rows:
        print(r)
