# Helix Guard — Railway Deployment Guide

This guide explains how to deploy **Helix Guard** on [Railway](https://railway.app/).

---

## Why Railway?
Railway is extremely fast and reliable for full-stack Python + Node apps:
- **Instant Deployment**: Deploys via Nixpacks build engine with zero configuration.
- **Single Service Full-Stack**: Helix Guard automatically bundles and serves both the React frontend and FastAPI backend on a single Railway domain!
- **Auto-Healing**: Automatic restarts on failure and zero-downtime rolling deploys.

---

## Method 1: Deploy from GitHub (Easiest & Recommended)

### Step 1: Push Code to GitHub
Ensure your latest code is committed and pushed to your GitHub repository:
```bash
git add .
git commit -m "feat: configure Railway deployment"
git push origin main
```

### Step 2: Create Railway Project
1. Go to [railway.app](https://railway.app/) and log in (or sign up with GitHub).
2. Click **+ New Project** in your dashboard.
3. Select **Deploy from GitHub repo**.
4. Choose your `helix-guard` (or `genome-firewall`) repository.

### Step 3: Railway Auto-Detection
Railway will automatically:
- Read [`railway.json`](./railway.json) and [`nixpacks.toml`](./nixpacks.toml).
- Detect Python 3.12 and install dependencies from `requirements.txt`.
- Start the server using:
  ```bash
  uvicorn app.main:app --host 0.0.0.0 --port $PORT
  ```
- Build and serve the React UI and REST API together.

### Step 4: Generate Your Public Domain
1. In your Railway project canvas, click on your service card.
2. Go to the **Settings** tab.
3. Scroll down to the **Networking** section.
4. Under **Public Networking**, click **Generate Domain**.
5. Railway will give you a live HTTPS domain, for example:
   ```
   https://helix-guard-production.up.railway.app
   ```
6. Open this URL in your browser:
   - `/` → Loads the **Helix Guard** React frontend with rotating 3D DNA canvas!
   - `/api/health` → Returns `{"status": "UP"}`
   - `/docs` → Interactive OpenAPI / Swagger documentation

---

## Method 2: Deploy Using Railway CLI

If you prefer deploying directly from your computer without GitHub:

1. Install Railway CLI:
   ```bash
   npm install -g @railway/cli
   ```
2. Login to your Railway account:
   ```bash
   railway login
   ```
3. Initialize a project inside the repository folder:
   ```bash
   railway init
   ```
4. Deploy the code:
   ```bash
   railway up
   ```
5. Generate a public domain:
   ```bash
   railway domain
   ```

---

## Method 3: Split Deployment (Backend on Railway, Frontend on Netlify)

If you prefer hosting the React frontend on Netlify's global CDN:

1. Deploy the backend on Railway following **Method 1**.
2. Copy your public Railway URL (e.g. `https://helix-guard-backend.up.railway.app`).
3. Deploy the frontend to Netlify:
   - Drag & drop `helix-guard-frontend-dist.zip` onto [app.netlify.com/drop](https://app.netlify.com/drop).
   - In Netlify **Site configuration → Environment variables**, add:
     - `VITE_API_URL` = `https://helix-guard-backend.up.railway.app`
4. The Netlify frontend will communicate securely with your Railway backend API over HTTPS.

---

## Environment Variables on Railway (Optional)

In Railway's **Variables** tab, you can customize these optional settings:
- `MAX_FILE_SIZE_MB`: Maximum FASTA file upload size (default: `50`).
- `PORT`: Automatically set by Railway (default: injected dynamically).
- `PYTHONUNBUFFERED`: `1` (for instant realtime log streaming).
