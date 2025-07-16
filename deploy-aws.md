# AWS Deployment Guide for Lipid Rendering

## 🚀 Option 1: AWS Amplify (Recommended)

### Prerequisites
- AWS Account
- GitHub repository (✅ Already set up)

### Steps:

1. **Go to AWS Amplify Console**
   - Navigate to https://console.aws.amazon.com/amplify/
   - Click "Get Started" under "Deploy"

2. **Connect Repository**
   - Choose "GitHub"
   - Select your repository: `Mahmoudimahyar/lipid-rendering`
   - Branch: `main`

3. **Configure Build Settings**
   - Amplify will auto-detect the `amplify.yml` file
   - Root directory: leave empty (it will use the amplify.yml configuration)

4. **Review and Deploy**
   - Click "Save and Deploy"
   - Your app will be available at: `https://[app-id].amplifyapp.com`

### Environment Variables (Optional)
In Amplify Console → App Settings → Environment Variables:
```
VITE_APP_NAME=Lipid Rendering
VITE_APP_VERSION=1.0.0
```

---

## ⚡ Option 2: Amazon S3 + CloudFront

### Prerequisites
- AWS CLI installed
- AWS credentials configured

### Steps:

1. **Build the Application**
   ```bash
   cd lipid_viewer
   npm run build:aws
   ```

2. **Create S3 Bucket**
   ```bash
   aws s3 mb s3://lipid-rendering-app --region us-east-1
   ```

3. **Upload Files**
   ```bash
   aws s3 sync dist/ s3://lipid-rendering-app --delete
   ```

4. **Configure S3 for Static Website**
   ```bash
   aws s3 website s3://lipid-rendering-app --index-document index.html --error-document index.html
   ```

5. **Set Bucket Policy for Public Access**
   ```json
   {
     "Version": "2012-10-17",
     "Statement": [
       {
         "Effect": "Allow",
         "Principal": "*",
         "Action": "s3:GetObject",
         "Resource": "arn:aws:s3:::lipid-rendering-app/*"
       }
     ]
   }
   ```

6. **Create CloudFront Distribution** (Optional but recommended)
   - Origin: Your S3 bucket website endpoint
   - Default behavior: Redirect HTTP to HTTPS
   - Custom error pages: 404 → /index.html (for SPA routing)

---

## 🐳 Option 3: AWS Elastic Beanstalk

### Prerequisites
- EB CLI installed
- Dockerfile created

### Steps:

1. **Initialize Elastic Beanstalk**
   ```bash
   cd lipid_viewer
   eb init lipid-rendering --platform node.js --region us-east-1
   ```

2. **Create Environment**
   ```bash
   eb create production
   ```

3. **Deploy**
   ```bash
   npm run build:aws
   eb deploy
   ```

---

## 📊 Cost Estimation

| Service | Monthly Cost (Est.) | Best For |
|---------|-------------------|----------|
| **Amplify** | $1-5 | Simplicity + CI/CD |
| **S3 + CloudFront** | $0.50-2 | Cost optimization |
| **Elastic Beanstalk** | $10-20 | Complex applications |

---

## 🔧 Production Optimizations

The following optimizations are already configured:

1. **Bundle Splitting**: Molecular libraries separated into chunks
2. **Asset Optimization**: Images and files optimized
3. **Caching**: Proper cache headers
4. **SPA Routing**: Redirects configured for React Router
5. **Security Headers**: CSP and CORS configured

---

## 🚨 Important Notes

1. **Molecular Libraries**: Your app uses large 3D molecular visualization libraries
   - Consider enabling gzip compression
   - Monitor bundle sizes (current molecular chunk ~2-3MB)

2. **Memory Usage**: 3D rendering can be memory-intensive
   - Consider Lambda@Edge for mobile optimization
   - Add loading states for large molecules

3. **CORS**: If using external molecular data APIs
   - Configure proper CORS headers
   - Use API Gateway if needed

---

## 🔍 Monitoring & Analytics

Add these to your deployment:

1. **CloudWatch**: For performance monitoring
2. **AWS X-Ray**: For request tracing
3. **AWS Pinpoint**: For user analytics
4. **Real User Monitoring**: For 3D performance tracking 