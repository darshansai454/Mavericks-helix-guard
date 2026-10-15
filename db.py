import sqlite3
import json
import os
from typing import List, Dict, Any, Optional
from datetime import datetime

DB_PATH = os.path.join(os.path.dirname(__file__), "genome_firewall.db")

def get_db():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    conn = get_db()
    cursor = conn.cursor()
    
    # Analyses table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS analyses (
        id TEXT PRIMARY KEY,
        created_at TEXT NOT NULL,
        organism TEXT NOT NULL,
        filename TEXT NOT NULL,
        file_sha256 TEXT NOT NULL,
        status TEXT NOT NULL,
        pipeline_version TEXT NOT NULL,
        model_version TEXT NOT NULL,
        summary_allow INTEGER DEFAULT 0,
        summary_block INTEGER DEFAULT 0,
        summary_review INTEGER DEFAULT 0,
        summary_total INTEGER DEFAULT 0,
        full_json TEXT NOT NULL
    );
    """)

    # Audit events table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS audit_events (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        analysis_id TEXT NOT NULL,
        timestamp TEXT NOT NULL,
        stage TEXT NOT NULL,
        event TEXT NOT NULL,
        details TEXT NOT NULL,
        FOREIGN KEY (analysis_id) REFERENCES analyses (id) ON DELETE CASCADE
    );
    """)

    # Settings table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS settings (
        key TEXT PRIMARY KEY,
        value TEXT NOT NULL
    );
    """)

    conn.commit()
    conn.close()

# Auto-initialize database on import
init_db()

def save_analysis(analysis_data: Dict[str, Any]):
    conn = get_db()
    cursor = conn.cursor()
    
    summary = analysis_data.get("summary_counts", {})
    cursor.execute("""
        INSERT OR REPLACE INTO analyses (
            id, created_at, organism, filename, file_sha256, status,
            pipeline_version, model_version, summary_allow, summary_block,
            summary_review, summary_total, full_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        analysis_data["analysis_id"],
        analysis_data["created_at"],
        analysis_data["organism"],
        analysis_data["filename"],
        analysis_data["file_sha256"],
        analysis_data["status"],
        analysis_data["pipeline_version"],
        analysis_data["model_version"],
        summary.get("allow", 0),
        summary.get("block", 0),
        summary.get("review", 0),
        summary.get("total", 0),
        json.dumps(analysis_data)
    ))

    # Save audit trail
    for event in analysis_data.get("audit_trail", []):
        cursor.execute("""
            INSERT INTO audit_events (analysis_id, timestamp, stage, event, details)
            VALUES (?, ?, ?, ?, ?)
        """, (
            analysis_data["analysis_id"],
            event["timestamp"],
            event["stage"],
            event["event"],
            event["details"]
        ))

    conn.commit()
    conn.close()

def get_analysis_by_id(analysis_id: str) -> Optional[Dict[str, Any]]:
    conn = get_db()
    cursor = conn.cursor()
    cursor.execute("SELECT full_json FROM analyses WHERE id = ?", (analysis_id,))
    row = cursor.fetchone()
    conn.close()
    if row:
        return json.loads(row["full_json"])
    return None

def list_analyses() -> List[Dict[str, Any]]:
    conn = get_db()
    cursor = conn.cursor()
    cursor.execute("""
        SELECT id as analysis_id, created_at, organism, filename, file_sha256, status,
               summary_allow, summary_block, summary_review, summary_total
        FROM analyses
        ORDER BY created_at DESC
    """)
    rows = cursor.fetchall()
    conn.close()
    
    results = []
    for r in rows:
        results.append({
            "analysis_id": r["analysis_id"],
            "created_at": r["created_at"],
            "organism": r["organism"],
            "filename": r["filename"],
            "file_sha256": r["file_sha256"],
            "status": r["status"],
            "verdicts": {
                "allow": r["summary_allow"],
                "block": r["summary_block"],
                "review": r["summary_review"],
                "total": r["summary_total"]
            }
        })
    return results

def delete_analysis_by_id(analysis_id: str) -> bool:
    conn = get_db()
    cursor = conn.cursor()
    cursor.execute("DELETE FROM audit_events WHERE analysis_id = ?", (analysis_id,))
    cursor.execute("DELETE FROM analyses WHERE id = ?", (analysis_id,))
    deleted = cursor.rowcount > 0
    conn.commit()
    conn.close()
    return deleted
