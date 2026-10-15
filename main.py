import os
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager

from .database.db import init_db
from .api.endpoints import router as api_router

@asynccontextmanager
async def lifespan(app: FastAPI):
    # Initialize SQLite database on startup
    init_db()
    yield

app = FastAPI(
    title="Helix Guard — AMR Prediction System",
    description=(
        "An Explainable Antimicrobial Resistance Prediction & Decision-Support System. "
        "RESEARCH-USE ONLY: Predictions are not clinical recommendations and must not be used alone "
        "to select or withhold antibiotic treatment."
    ),
    version="2.4.0",
    lifespan=lifespan
)

# Enable CORS for local Vite development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount API endpoints
app.include_router(api_router, prefix="/api")

# Check for pre-built frontend distribution
candidates_dist = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "frontend", "dist")),
    os.path.abspath(os.path.join(os.getcwd(), "frontend", "dist")),
    os.path.abspath(os.path.join(os.getcwd(), "dist")),
]
frontend_dist = next((d for d in candidates_dist if os.path.exists(os.path.join(d, "index.html"))), None)

if frontend_dist:
    from fastapi.staticfiles import StaticFiles
    app.mount("/", StaticFiles(directory=frontend_dist, html=True), name="frontend")
else:
    @app.get("/")
    def root():
        return {
            "system": "Helix Guard",
            "version": "2.4.0",
            "status": "OPERATIONAL",
            "documentation": "/docs",
            "disclaimer": (
                "Helix Guard is a research prototype for antimicrobial resistance prediction. "
                "Its predictions are not clinical recommendations."
            )
        }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app.main:app", host="127.0.0.1", port=8000, reload=True)
