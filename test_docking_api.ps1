# Test docking API with a simple request

$headers = @{
    'Content-Type' = 'application/json'
}

$body = @{
    ligand_smiles = "CCO"
    receptor_pdb_id = "1CRN"
    center_x = 0.0
    center_y = 0.0
    center_z = 0.0
    size_x = 20.0
    size_y = 20.0
    size_z = 20.0
    exhaustiveness = 4
    num_modes = 3
    seed = 42
} | ConvertTo-Json

Write-Host "Testing docking API with simple request..." -ForegroundColor Cyan
Write-Host "Request body: $body" -ForegroundColor Yellow

try {
    $response = Invoke-WebRequest -Uri "http://localhost:8000/api/dock/run" -Method POST -Body $body -Headers $headers -UseBasicParsing
    Write-Host "✅ Success! Status: $($response.StatusCode)" -ForegroundColor Green
    Write-Host "Response: $($response.Content)" -ForegroundColor White
} catch {
    Write-Host "❌ Error: $($_.Exception.Message)" -ForegroundColor Red
    if ($_.Exception.Response) {
        $errorContent = $_.Exception.Response.Content.ReadAsStringAsync().Result
        Write-Host "Error details: $errorContent" -ForegroundColor Yellow
    }
}
