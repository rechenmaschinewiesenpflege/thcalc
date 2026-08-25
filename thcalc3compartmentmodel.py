import math

def simulate_thc_clearance(
    gender, age, weight, height, body_fat, activity, weight_trend, hydration, diet,
    years_used, days_per_week, grams_per_day, last_break, product_source,
    method, tobacco_mix, alcohol_mix, thc_percent, cbd_level, dose_amount,
    drag_size, hold_time, session_duration,
    high_duration, edible_sens, grogginess, fasting_state,
    meds_inhibitors, post_activity
):
    # --- 1. CLINICALLY SCALED BIOMETRICS ---
    if body_fat <= 0:
        bmi = weight / ((height / 100) ** 2)
        gender_factor = 1 if gender == "m" else 0
        body_fat = max(5.0, (1.20 * bmi) + (0.23 * age) - (10.8 * gender_factor) - 5.4)
    
    # Apparent Initial Distribution Volume (V1_apparent)
    v1_apparent = weight * 0.80 * {1: 1.05, 2: 1.0, 3: 0.9}.get(hydration, 1.0) 

    # --- 2. ENZYME PHENOTYPING (CYP2C9 / CYP3A4 Proxies) ---
    cyp_score = 3.0
    cyp_score += {1: -1.0, 2: -0.5, 3: 0, 4: 0.5, 5: 1.0}.get(high_duration, 0)
    cyp_score += {1: 1.0, 2: 0.5, 3: 0, 4: -0.5, 5: -1.0}.get(edible_sens, 0)
    cyp_score += {1: -1.0, 2: -0.5, 3: 0, 4: 0.5, 5: 1.0}.get(grogginess, 0)
    cyp_phenotype = max(1.0, min(5.0, cyp_score))

    # --- 3. BIOAVAILABILITY (F) & ABSORPTION ---
    f = 0.0
    ka = 12.0
    if method == 1: f = 0.25      
    elif method == 2: f = 0.30    
    elif method == 3: f = 0.40    
    elif method == 4: f = 0.45    
    elif method == 5: f = 0.10    
    elif method == 6: f = 0.20    
    
    if method in [1,2,3,4]:
        f *= {1: 0.7, 2: 1.0, 3: 1.3}.get(drag_size, 1.0)
        f *= {1: 0.8, 2: 1.0, 3: 1.2, 4: 1.4}.get(hold_time, 1.0)
        if tobacco_mix == 1: f *= 0.85 
        if alcohol_mix == 1: f *= 1.20 
    else:
        if fasting_state == 1: ka = 1.5; f *= 0.8 
        elif fasting_state == 2: ka = 1.0; f *= 1.0
        elif fasting_state == 3: ka = 0.6; f *= 1.4 
        
    actual_thc_percent = thc_percent * (0.85 if product_source == 2 else 1.0)
    dose_mg = dose_amount * 1000 * (actual_thc_percent / 100.0)
    absorbed_mcg = dose_mg * f * 1000

    # --- 4. CALIBRATED TRI-EXPONENTIAL CONSTANTS ---
    alpha = math.log(2) / 0.15 
    
    # Beta: Hepatic clearance
    base_t_half_beta = 3.0 - (cyp_phenotype - 3.0) * 0.4 
    
    # CBD is a competitive inhibitor of CYP enzymes
    if cbd_level == 1: base_t_half_beta *= 1.3 
    elif cbd_level == 2: base_t_half_beta *= 1.15
    
    # Detailed Medication Interactions (CYP2C9 / CYP3A4)
    if meds_inhibitors == 2: 
        base_t_half_beta *= 1.55 # Strong Inhibitors massively delay clearance
    elif meds_inhibitors == 3: 
        base_t_half_beta *= 0.60 # Strong Inducers accelerate clearance
    elif meds_inhibitors == 4: 
        base_t_half_beta *= 1.25 # Mild Inhibitors gently delay clearance
        
    beta = math.log(2) / base_t_half_beta
    
    # Gamma: Terminal lipid release half-life
    terminal_days = 3.0 + (body_fat / 10.0) + min(4.0, (years_used * days_per_week) / 10.0)
    gamma = math.log(2) / (terminal_days * 24.0)

    # Re-balanced Coefficients (90% distributes rapidly)
    C0 = absorbed_mcg / v1_apparent
    A = C0 * 0.90      
    B = C0 * 0.085     
    C_acute = C0 * 0.015 

    # --- 5. CHRONIC LIPID BASELINE (C_chronic) ---
    chronicity_factor = min(1.0, years_used / 3.0) * (days_per_week / 7.0)
    daily_mg = grams_per_day * 1000 * (actual_thc_percent / 100.0) * f
    break_reduction = max(0.0, 1.0 - (last_break / 30.0))
    
    lipid_storage_mg = daily_mg * chronicity_factor * break_reduction * (body_fat / 15.0)
    baseline_c0 = lipid_storage_mg / (weight * 0.5) 
    
    lipolysis_mod = 1.0
    if weight_trend == 1: lipolysis_mod *= 1.3 
    elif weight_trend == 3: lipolysis_mod *= 0.8 
    if post_activity == 3: lipolysis_mod *= 1.5 
    elif post_activity == 1: lipolysis_mod *= 0.9
    
    baseline_c0 *= lipolysis_mod
    C_terminal_total = C_acute + baseline_c0

    # --- 6. TIMELINE SIMULATION (30 DAYS) ---
    simulation_hours = 30 * 24 
    timeline = []
    threshold_hour = -1.0
    
    rsd = 0.25 + (0.10 if product_source == 2 else 0)

    for step in range(1, (simulation_hours * 2) + 1): 
        t = step / 2.0
        
        absorption_factor = (ka / (ka - alpha)) if ka != alpha else 1.0
        
        conc = absorption_factor * (
            A * math.exp(-alpha * t) + 
            B * math.exp(-beta * t) + 
            C_terminal_total * math.exp(-gamma * t)
        ) - (C0 * math.exp(-ka * t))
        
        conc = max(0.0, conc)
        
        ci_lower = max(0.0, conc * (1.0 - 1.96 * rsd))
        ci_upper = conc * (1.0 + 1.96 * rsd)
        
        timeline.append({"hour": t, "conc": conc, "lower": ci_lower, "upper": ci_upper})
        
        if threshold_hour == -1.0 and conc < 3.5:
            threshold_hour = t

    return timeline, threshold_hour, rsd, cyp_phenotype, f, C0, baseline_c0

if __name__ == "__main__":
    print("\n" + "="*80)
    print(" 🌿 CLINICALLY CALIBRATED 3-COMPARTMENT THC MODEL (3.5 ng/ml)")
    print("="*80)
    print(" Tip: Type 'b' to go back and correct mistakes.\n")

    prompts = [
        # --- BIOMETRICS ---
        {"key": "gender", "msg": "1. Biological Gender (m/f/x)", "type": "str", "default": "m"},
        {"key": "age", "msg": "2. Age (years)", "type": "int", "default": 25},
        {"key": "weight", "msg": "3. Weight (kg)", "type": "float", "default": 75.0},
        {"key": "height", "msg": "4. Height (cm)", "type": "float", "default": 180.0},
        {"key": "body_fat", "msg": "5. Body Fat % (Enter 0 for auto-estimate)", "type": "float", "default": 0.0},
        
        # --- LIFESTYLE ---
        {"key": "activity", "msg": "6. General Activity Level (1: Sedentary, 3: Mod, 5: Athlete)", "type": "int", "default": 3, "min": 1, "max": 5},
        {"key": "weight_trend", "msg": "7. Recent weight trend? (1: Losing, 2: Stable, 3: Gaining)", "type": "int", "default": 2, "min": 1, "max": 3},
        {"key": "hydration", "msg": "8. Current Hydration (1: Dehydrated, 2: Normal, 3: Very Hydrated)", "type": "int", "default": 2, "min": 1, "max": 3},
        {"key": "diet", "msg": "9. General Diet (1: High Fat/Keto, 2: Balanced, 3: Low Fat)", "type": "int", "default": 2, "min": 1, "max": 3},
        
        # --- CHRONIC HISTORY ---
        {"key": "years", "msg": "10. Years of regular cannabis use?", "type": "float", "default": 3.0},
        {"key": "days_wk", "msg": "11. Days per week used on average?", "type": "int", "default": 4, "min": 0, "max": 7},
        {"key": "grams_day", "msg": "12. Grams per day (on usage days)?", "type": "float", "default": 0.5},
        {"key": "break", "msg": "13. Longest T-Break in the last 6 months (in days)?", "type": "int", "default": 0},
        {"key": "source", "msg": "14. Primary source (1: Medical/Dispensary, 2: Street/Homegrow)", "type": "int", "default": 1, "min": 1, "max": 2},
        
        # --- ACUTE SESSION DETAILS ---
        {"key": "thc_pct", "msg": "15. Estimated THC content (%)?", "type": "float", "default": 20.0},
        {"key": "cbd", "msg": "16. CBD content? (1: High 1:1, 2: Moderate, 3: Low/None)", "type": "int", "default": 3, "min": 1, "max": 3},
        {"key": "method", "msg": "17. Method (1: Joint, 2: Bong, 3: Dry Vape, 4: Cartridge, 5: Edible, 6: Sublingual)", "type": "int", "default": 1, "min": 1, "max": 6},
        {"key": "tobacco", "msg": "18. Mixed with tobacco/nicotine? (1: Yes, 2: No)", "type": "int", "default": 2, "min": 1, "max": 2},
        {"key": "alcohol", "msg": "19. Consumed alcohol alongside this session? (1: Yes, 2: No)", "type": "int", "default": 2, "min": 1, "max": 2},
        {"key": "dose", "msg": "20. Total amount consumed THIS session (grams/ml)?", "type": "float", "default": 0.3},
        
        # --- INHALATION MECHANICS ---
        {"key": "drag", "msg": "21. Average puff size? (1: Small/Sips, 2: Normal, 3: Deep lung)", "type": "int", "default": 2, "min": 1, "max": 3},
        {"key": "hold", "msg": "22. Smoke hold time? (1: Instant, 2: 1-2s, 3: 3-5s, 4: >5s)", "type": "int", "default": 2, "min": 1, "max": 4},
        {"key": "duration", "msg": "23. Session duration in minutes?", "type": "float", "default": 10.0},
        
        # --- ENZYME PROXIES (CYP) ---
        {"key": "high_dur", "msg": "24. Compared to friends, your high lasts: (1: Much longer -> 5: Much shorter)", "type": "int", "default": 3, "min": 1, "max": 5},
        {"key": "edible", "msg": "25. Sensitivity to edibles? (1: Too intense, 3: Normal, 5: Barely feel them)", "type": "int", "default": 3, "min": 1, "max": 5},
        {"key": "grogginess", "msg": "26. Next-day grogginess? (1: Always heavy, 3: Sometimes, 5: Never)", "type": "int", "default": 3, "min": 1, "max": 5},
        {"key": "fasting", "msg": "27. Stomach status? (1: Fasted, 2: Normal, 3: High Fat Meal)", "type": "int", "default": 2, "min": 1, "max": 3},
        
        # --- NEW DETAILED MEDICATION QUESTIONNAIRE ---
        {"key": "meds", "header": "\n--- 💊 Medication & Supplement Interactions (CYP2C9 & CYP3A4) ---\nLiver enzymes physically break down THC. Other drugs can block or speed this up:\n   [1] None of the below\n   [2] Strong Inhibitors (e.g., Fluconazole, Valproic Acid/Depakote, Ketoconazole, Amiodarone, Fluvoxamine)\n   [3] Strong Inducers (e.g., St. John's Wort, Carbamazepine, Rifampin, Phenobarbital)\n   [4] Mild/Moderate Inhibitors (e.g., Daily Grapefruit Juice, Omeprazole, Cimetidine)", 
         "msg": "28. Your choice (1-4)", "type": "int", "default": 1, "min": 1, "max": 4},
        
        {"key": "post_act", "msg": "29. Activity AFTER consuming? (1: Couch/Sleep, 2: Walking, 3: Heavy Gym)", "type": "int", "default": 1, "min": 1, "max": 3}
    ]

    answers = {}
    i = 0
    while i < len(prompts):
        p = prompts[i]
        
        if "header" in p: print(p["header"])
        user_input = input(f"{p['msg']} [Default: {p['default']}]: ").strip()

        if user_input.lower() == 'b':
            if i > 0:
                i -= 1
                print("\n" + "-"*40 + "\n ⏪ Going back...\n" + "-"*40)
            continue

        if user_input == "":
            answers[p["key"]] = p["default"]; i += 1; continue

        try:
            if p["type"] == "str": answers[p["key"]] = user_input.lower()
            elif p["type"] == "int":
                val = int(user_input)
                if "min" in p and (val < p["min"] or val > p["max"]): raise ValueError
                answers[p["key"]] = val
            elif p["type"] == "float": answers[p["key"]] = float(user_input.replace(",", "."))
            i += 1 
        except ValueError: print("❌ Invalid input. Try again.")

    timeline, threshold_hour, rsd, cyp_phenotype, f, c0, base_c0 = simulate_thc_clearance(
        answers["gender"], answers["age"], answers["weight"], answers["height"], answers["body_fat"], 
        answers["activity"], answers["weight_trend"], answers["hydration"], answers["diet"],
        answers["years"], answers["days_wk"], answers["grams_day"], answers["break"], answers["source"],
        answers["method"], answers["tobacco"], answers["alcohol"], answers["thc_pct"], answers["cbd"], 
        answers["dose"], answers["drag"], answers["hold"], answers["duration"],
        answers["high_dur"], answers["edible"], answers["grogginess"], answers["fasting"],
        answers["meds"], answers["post_act"]
    )

    print("\n" + "="*80)
    print(" 🔬 3-COMPARTMENT ELIMINATION ANALYSIS")
    print("="*80)
    print(f"Calculated Systemic Uptake (F): {round(f * 100, 1)}%")
    print(f"Hepatic Clearance Phenotype:  {cyp_phenotype}/5.0 (CYP2C9/3A4 Proxy Score)")
    print(f"Theoretical Initial Peak:     {round(c0, 1)} ng/ml (Pre-distribution)")
    print(f"Chronic Lipid Baseline:       {round(base_c0, 2)} ng/ml (Floor level)")
    
    print("\n" + "="*80)
    print(" 🚗 LEGAL DRIVING THRESHOLD (3.5 ng/ml)")
    print("="*80)
    if threshold_hour == -1.0:
        print("🚨 WARNING: Your chronic lipid baseline is so high that you will not drop below 3.5 ng/ml within 30 days of abstinence.")
    else:
        days = math.floor(threshold_hour / 24)
        hrs = threshold_hour % 24
        
        worst_case_hour = -1.0
        for data in timeline:
            if worst_case_hour == -1.0 and data["upper"] < 3.5:
                worst_case_hour = data["hour"]
                
        print(f"✅ Expected Mean Clearance Time:  {days} days and {hrs} hours")
        if worst_case_hour != -1.0:
            wc_days = math.floor(worst_case_hour / 24)
            wc_hrs = worst_case_hour % 24
            print(f"🛡️ Worst-Case (95% CI Upper):    {wc_days} days and {wc_hrs} hours")
        else:
            print("🛡️ Worst-Case (95% CI Upper):    Exceeds 30 days due to chronic lipid storage.")

    print("\n" + "="*80)
    print(" 📉 CLEARANCE TIMELINE (Blood Serum Concentration)")
    print("="*80)
    print(f"{'Time Passed':<15} | {'Mean (ng/ml)':<12} | {'95% CI Range (ng/ml)':<20} | {'Status'}")
    print("-" * 80)
    
    checkpoints = [1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0, 72.0, 168.0, 336.0]
    
    for data in timeline:
        t = data["hour"]
        if t in checkpoints or t == threshold_hour:
            conc = data["conc"]
            low = data["lower"]
            high = data["upper"]
            
            if t < 24: t_str = f"+ {t} hours"
            else: t_str = f"+ {t/24} days"
            
            if conc >= 3.5: status = "🔴 Impaired / Illegal"
            elif high >= 3.5: status = "🟡 Legal, but CI risks failing"
            elif conc >= 1.0: status = "🟢 Legal (Trace detectable)"
            else: status = "⚪ Completely Cleared"
            
            marker = " <== THRESHOLD CROSSED" if t == threshold_hour else ""
            
            print(f"{t_str:<15} | {round(conc, 2):<12} | {round(low, 2)} - {round(high, 2):<14} | {status}{marker}")
            
    print("-" * 80 + "\n")
