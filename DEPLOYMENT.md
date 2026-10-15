# Helix Guard — Deployment Guide (Netlify & Render)

This guide provides step-by-step instructions for deploying **Helix Guard**:
- **Backend API (FastAPI + ML Engine)**: Deployed on **Render** (Free Web Service)
- **Frontend App (React + Vite + Tailwind)**: Deployed on **Netlify** (Free Static Web App)

---

## Architecture Overview

```
                                      +---------------------------------------------+
                                      |          Netlify (Frontend CDN)             |
                                      |   https://helix-guard.netlify.app           |
                                      |   - React 18 + Vite + Tailwind CSS          |
                                      |   - 3D Animated DNA Laboratory Canvas       |
                                      |   - 18 Antibiotic Reference Catalog         |
                                      +----------------------+----------------------+
                                                             |
                                                             | HTTPS API Requests
                                                             | (VITE_API_URL)
                                                             v
                                      +---------------------------------------------+
                                      |          Render (Backend API)               |
                                      |   https://helix-guard-backend.onrender.com   |
                                      |   - FastAPI REST API                        |
                                      |   - Calibrated ML Prediction Engine         |
                                      |   - SQLite Audit & Report Storage           |
                                      +---------------------------------------------+
```

---

## Part 1: Deploy Backend on Render (Do this first)

You need your live Render backend URL before deploying the frontend so Netlify knows where to send API requests.

### Option A: 1-Click Render Blueprint (Recommended)
1. Push your project repository to GitHub or GitLab.
2. Log in to [Render Dashboard](https://dashboard.render.com/).
3. Click **New +** in the top navigation and select **Blueprint**.
4. Connect your GitHub repository.
5. Render will automatically read [`render.yaml`](./render.yaml) and configure the `helix-guard-backend` service.
6. Click **Apply**.
7. Wait 2-3 minutes for the build to finish. Your backend will be live at:
   ```
   https://helix-guard-backend-xxxx.onrender.com
   ```

---

### Option B: Manual Web Service Setup on Render
If you prefer configuring manually through the Render UI:

1. In [Render Dashboard](https://dashboard.render.com/), click **New +** -> **Web Service**.
2. Connect your GitHub repository.
3. Configure the following fields:
   - **Name**: `helix-guard-backend`
   - **Region**: Any close to you (e.g., *Oregon (US West)* or *Frankfurt (EU)*)
   - **Branch**: `main`
   - **Root Directory**: Leave blank (root of repository)
   - **Runtime**: `Python 3`
   - **Build Command**:
     ```bash
     pip install --upgrade pip && pip install -r backend/requirements.txt
     ```
   - **Start Command**:
     ```bash
     uvicorn backend.app.main:app --host 0.0.0.0 --port $PORT
     ```
   - **Instance Type**: `Free`
4. Expand **Advanced** -> **Add Environment Variable**:
   - `PYTHON_VERSION`: `3.12.0`
5. Click **Create Web Service**.
6. When deployment finishes, copy your service URL (e.g., `https://helix-guard-backend.onrender.com`).
7. **Verify it is working**: Open `https://YOUR-RENDER-URL.onrender.com/api/health` in your browser. You should see:
   ```json
   {
     "status": "UP",
     "service": "Helix Guard Prediction Engine",
     "version": "2.4.0",
     "database_connected": true
   }
   ```

> **Note on Render Free Tier**: On free tier, Render puts web services to sleep after 15 minutes of inactivity. When a request arrives, it may take 30-50 seconds to wake up (cold start).

---

## Part 2: Deploy Frontend on Netlify

Now deploy your React frontend and link it to your live Render backend.

### Step 1: Connect Repository to Netlify
1. Log in to your [Netlify Dashboard](https://app.netlify.com/).
2. Click **Add new site** -> **Import an existing project**.
3. Choose **GitHub** and authorize access to your repository.

### Step 2: Configure Build Settings
Netlify will automatically detect the settings from [`netlify.toml`](./netlify.toml), but verify the following values:
- **Base directory**: `frontend`
- **Build command**: `npm run build`
- **Publish directory**: `frontend/dist` (or `dist` if base directory is set to `frontend`)

### Step 3: Add the Backend API URL
1. Still on the configuration page, click **Add environment variables**.
2. Add the variable:
   - **Key**: `VITE_API_URL`
   - **Value**: `https://YOUR-RENDER-BACKEND.onrender.com` *(paste your actual Render backend URL without trailing slash)*
3. Click **Deploy Helix Guard**.

### Step 4: Verification
1. Once Netlify completes the build (usually ~45 seconds), click your new Netlify URL (e.g., `https://helix-guard.netlify.app`).
2. Verify:
   - The laboratory hero with rotating 3D DNA canvas loads smoothly.
   - Click **Antibiotic Reference** in the navigation bar to see all 18 clinical compounds.
   - Click **Upload Genome** -> choose the benchmark *MDR Klebsiella (NDM-1)* sample -> click **Run Helix Guard Analysis**.
   - Verify the light-themed report appears with the **"GENOMIC PREDICTION (NOT EXPERIMENTAL AST)"** badge.

---

## Part 3: SPA Routing & CORS Support

Everything is pre-configured out-of-the-box:
1. **Single Page Application Routing**: Both [`netlify.toml`](./netlify.toml) and [`frontend/public/_redirects`](./frontend/public/_redirects) contain the rewrite `/* -> /index.html 200` to prevent 404 errors when refreshing on `/dashboard`, `/methodology`, etc.
2. **CORS (Cross-Origin Resource Sharing)**: The FastAPI backend in [`backend/app/main.py`](./backend/app/main.py) allows requests from any origin (`allow_origins=["*"]`), ensuring Netlify can securely fetch API data without browser CORS errors.

---

## Part 4: Local Development Quick Start

To run both services locally on your machine before or alongside deployment:

### Backend:
```bash
cd backend
python -m venv venv
# Windows:
.\venv\Scripts\activate
# Linux/macOS:
source venv/bin/activate

pip install -r requirements.txt
uvicorn backend.app.main:app --host 127.0.0.1 --port 8000 --reload
```

### Frontend:
```bash
cd frontend
npm install
npm run dev
```
Open [http://127.0.0.1:5173](http://127.0.0.1:5173) in your browser.
