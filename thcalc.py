import math
import matplotlib.pyplot as plt

def simulate_thc_clearance(
    gender, age, weight, height, body_fat, activity, weight_trend, hydration, diet,
    years_used, days_per_week, grams_per_day, last_break, product_source,
    method, tobacco_mix, alcohol_mix, thc_percent, cbd_level, dose_amount,
    drag_size, hold_time, session_duration,
    high_duration, edible_sens, grogginess, fasting_state,
    med_list, post_activity
):
    # --- 1. CLINICALLY SCALED BIOMETRICS ---
    if body_fat <= 0:
        bmi = weight / ((height / 100) ** 2)
        gender_factor = 1 if gender == "m" else 0
        body_fat = max(5.0, (1.20 * bmi) + (0.23 * age) - (10.8 * gender_factor) - 5.4)
    
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
    
    if method in [1, 2, 3, 4]:
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
    base_t_half_beta = 3.0 - (cyp_phenotype - 3.0) * 0.4 
    
    if cbd_level == 1: base_t_half_beta *= 1.3 
    elif cbd_level == 2: base_t_half_beta *= 1.15
    
    med_warnings = []
    med_factor = 1.0
    
    if any(m in med_list for m in [1, 2, 3, 4]):
        med_factor *= 1.60
        med_warnings.append("Strong CYP2C9 Inhibitor detected. Clearance heavily delayed.")
    if any(m in med_list for m in [5, 6, 7]):
        med_factor *= 1.25
        med_warnings.append("Moderate CYP2C9 Inhibitor detected. Clearance moderately delayed.")
    if any(m in med_list for m in [8, 9, 10]):
        med_factor *= 1.40
        med_warnings.append("Strong CYP3A4 Inhibitor detected. Secondary clearance pathway blocked.")
    if any(m in med_list for m in [11, 12, 13, 14]):
        med_factor *= 0.55
        med_warnings.append("Strong Hepatic Inducer detected. Clearance heavily accelerated.")

    base_t_half_beta *= med_factor
    beta = math.log(2) / base_t_half_beta
    
    terminal_days = 3.0 + (body_fat / 10.0) + min(4.0, (years_used * days_per_week) / 10.0)
    gamma = math.log(2) / (terminal_days * 24.0)

    C0 = absorbed_mcg / v1_apparent
    A = C0 * 0.90      
    B = C0 * 0.085     
    C_acute = C0 * 0.015 

    # --- 5. CHRONIC LIPID BASELINE ---
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

    # --- 6. TIMELINE SIMULATION ---
    simulation_hours = 30 * 24 
    timeline = {}
    threshold_hour = -1.0
    cleared_hour = -1.0
    
    rsd = 0.25 + (0.10 if product_source == 2 else 0)

    for step in range(1, (simulation_hours * 2) + 1): 
        t = round(step / 2.0, 1)
        absorption_factor = (ka / (ka - alpha)) if ka != alpha else 1.0
        
        conc = absorption_factor * (
            A * math.exp(-alpha * t) + 
            B * math.exp(-beta * t) + 
            C_terminal_total * math.exp(-gamma * t)
        ) - (C0 * math.exp(-ka * t))
        
        conc = max(0.0, conc)
        ci_lower = max(0.0, conc * (1.0 - 1.96 * rsd))
        ci_upper = conc * (1.0 + 1.96 * rsd)
        
        timeline[t] = {"conc": conc, "lower": ci_lower, "upper": ci_upper}
        
        if threshold_hour == -1.0 and conc < 3.5:
            threshold_hour = t
        if cleared_hour == -1.0 and conc < 1.0:
            cleared_hour = t

    return timeline, threshold_hour, cleared_hour, rsd, cyp_phenotype, f, C0, baseline_c0, med_warnings

def run_medication_questionnaire():
    med_catalog = {
        1: "Fluconazole / Diflucan (Strong CYP2C9 Inhibitor)",
        2: "Valproic Acid / Depakote (Strong CYP2C9 Inhibitor)",
        3: "Amiodarone / Cordarone (Strong CYP2C9 Inhibitor)",
        4: "Miconazole / Voriconazole (Strong CYP2C9 Inhibitors)",
        5: "Fluvoxamine / Luvox (Moderate CYP2C9 Inhibitor)",
        6: "Sertraline / Zoloft (Mild/Moderate CYP2C9 Inhibitor)",
        7: "Omeprazole / Pantoprazole / Cimetidine (Mild CYP2C9 Inhibitor)",
        8: "Ketoconazole / Itraconazole (Strong CYP3A4 Inhibitors)",
        9: "Clarithromycin / Erythromycin (Strong CYP3A4 Inhibitors)",
        10: "Grapefruit Juice Extract (CYP3A4 Inhibitor)",
        11: "St. John's Wort / Johanniskraut (Strong Inducer)",
        12: "Rifampin / Rifampicin (Potent Inducer)",
        13: "Carbamazepine / Tegretol (Strong Inducer)",
        14: "Phenytoin / Phenobarbital (Strong Inducers)"
    }
    selected_meds = set()
    print("\n" + "="*80 + "\n 💊 MEDICATION & SUBSTANCE INTERACTION QUESTIONNAIRE\n" + "="*80)
    for num, name in med_catalog.items(): print(f"   [{num:2d}] {name}")
    print("   [ 0] Done / Finish selection\n" + "-" * 80)
    
    while True:
        current_selection_str = ", ".join(str(m) for m in sorted(selected_meds)) if selected_meds else "None"
        print(f"\nCurrent Selection: [{current_selection_str}]")
        raw_val = input("Enter medication number to add/remove (or 0 to finish): ").strip()
        if raw_val == "0" or raw_val == "": break
        try:
            choice = int(raw_val)
            if choice in med_catalog:
                if choice in selected_meds:
                    selected_meds.remove(choice)
                    print(f"➖ Removed: {med_catalog[choice]}")
                else:
                    selected_meds.add(choice)
                    print(f"➕ Added: {med_catalog[choice]}")
            else: print("❌ Invalid number.")
        except ValueError: print("❌ Invalid input.")
    return list(selected_meds)


def plot_pharmacokinetics(timeline, threshold_hour, cleared_hour, baseline_c0):
    """Generates a dynamic matplotlib chart focused on the legal threshold crossing."""
    
    # 1. Determine dynamic X-axis limit
    if threshold_hour != -1.0:
        max_t = threshold_hour * 1.5
        if cleared_hour != -1.0:
            max_t = max(max_t, cleared_hour * 1.1)
    else:
        max_t = 72.0 # Default fallback if 3.5 limit is never crossed
        
    # Extract data for plotting within the dynamic range
    x_hours = []
    y_mean = []
    y_lower = []
    y_upper = []
    
    for t, data in sorted(timeline.items()):
        if t <= max_t:
            x_hours.append(t)
            y_mean.append(data["conc"])
            y_lower.append(data["lower"])
            y_upper.append(data["upper"])

    plt.figure(figsize=(12, 7))
    
    # Plot Mean and 95% Confidence Interval
    plt.plot(x_hours, y_mean, color='#1f77b4', linewidth=2.5, label='Mean Plasma THC (ng/ml)')
    plt.fill_between(x_hours, y_lower, y_upper, color='#1f77b4', alpha=0.2, label='95% Confidence Interval')

    # Add Legal Threshold Lines
    plt.axhline(y=3.5, color='red', linestyle='--', linewidth=2, label='Legal Driving Threshold (3.5 ng/ml)')
    plt.axhline(y=1.0, color='green', linestyle=':', linewidth=2, label='Fully Cleared (< 1.0 ng/ml)')

    # Add Exact Crossing Annotations
    if threshold_hour != -1.0:
        plt.plot(threshold_hour, 3.5, 'ro', markersize=8)
        plt.annotate(f'3.5 ng/ml crossed at\n{threshold_hour:.1f} hours', 
                     xy=(threshold_hour, 3.5), xytext=(threshold_hour + (max_t*0.05), 4.5),
                     arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=6),
                     fontsize=10, fontweight='bold', color='red')
                     
    if cleared_hour != -1.0 and cleared_hour <= max_t:
        plt.plot(cleared_hour, 1.0, 'go', markersize=8)
        plt.annotate(f'1.0 ng/ml crossed at\n{cleared_hour:.1f} hours', 
                     xy=(cleared_hour, 1.0), xytext=(cleared_hour + (max_t*0.05), 1.5),
                     arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=6),
                     fontsize=10, fontweight='bold', color='green')

    # Dynamic Y-Axis Capping (Prevents initial peak from compressing the chart)
    y_max = max(15.0, baseline_c0 * 2.5) 
    plt.ylim(0, y_max)
    plt.xlim(0, max_t)

    # Styling & Labels
    plt.title("3-Compartment THC Pharmacokinetics & Legal Threshold Tracking", fontsize=14, fontweight='bold', pad=15)
    plt.xlabel("Time Since Consumption (Hours)", fontsize=12, fontweight='bold')
    plt.ylabel("Active Plasma THC (ng/ml)", fontsize=12, fontweight='bold')
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.legend(loc='upper right', fontsize=10, framealpha=0.9)
    plt.tight_layout()
    
    # Display Chart
    plt.show()


if __name__ == "__main__":
    print("\n" + "="*80)
    print(" 🌿 3-COMPARTMENT THC MODEL WITH MATPLOTLIB VISUALIZATION")
    print("="*80)
    print(" Tip: Type 'b' to go back at any prompt.\n")

    prompts = [
        {"key": "gender", "msg": "1. Biological Gender (m/f/x)", "type": "str", "default": "m"},
        {"key": "age", "msg": "2. Age (years)", "type": "int", "default": 25},
        {"key": "weight", "msg": "3. Weight (kg)", "type": "float", "default": 75.0},
        {"key": "height", "msg": "4. Height (cm)", "type": "float", "default": 180.0},
        {"key": "body_fat", "msg": "5. Body Fat % (Enter 0 for auto-estimate)", "type": "float", "default": 0.0},
        {"key": "activity", "msg": "6. General Activity Level (1: Sedentary, 3: Mod, 5: Athlete)", "type": "int", "default": 3, "min": 1, "max": 5},
        {"key": "weight_trend", "msg": "7. Recent weight trend? (1: Losing, 2: Stable, 3: Gaining)", "type": "int", "default": 2, "min": 1, "max": 3},
        {"key": "hydration", "msg": "8. Current Hydration (1: Dehydrated, 2: Normal, 3: Very Hydrated)", "type": "int", "default": 2, "min": 1, "max": 3},
        {"key": "diet", "msg": "9. General Diet (1: High Fat/Keto, 2: Balanced, 3: Low Fat)", "type": "int", "default": 2, "min": 1, "max": 3},
        {"key": "years", "msg": "10. Years of regular cannabis use?", "type": "float", "default": 3.0},
        {"key": "days_wk", "msg": "11. Days per week used on average?", "type": "int", "default": 4, "min": 0, "max": 7},
        {"key": "grams_day", "msg": "12. Grams per day (on usage days)?", "type": "float", "default": 0.5},
        {"key": "break", "msg": "13. Longest T-Break in the last 6 months (in days)?", "type": "int", "default": 0},
        {"key": "source", "msg": "14. Primary source (1: Medical/Dispensary, 2: Street/Homegrow)", "type": "int", "default": 1, "min": 1, "max": 2},
        {"key": "thc_pct", "msg": "15. Estimated THC content (%)?", "type": "float", "default": 20.0},
        {"key": "cbd", "msg": "16. CBD content? (1: High 1:1, 2: Moderate, 3: Low/None)", "type": "int", "default": 3, "min": 1, "max": 3},
        {"key": "method", "msg": "17. Method (1: Joint, 2: Bong, 3: Dry Vape, 4: Cartridge, 5: Edible, 6: Sublingual)", "type": "int", "default": 1, "min": 1, "max": 6},
        {"key": "tobacco", "msg": "18. Mixed with tobacco/nicotine? (1: Yes, 2: No)", "type": "int", "default": 2, "min": 1, "max": 2},
        {"key": "alcohol", "msg": "19. Consumed alcohol alongside this session? (1: Yes, 2: No)", "type": "int", "default": 2, "min": 1, "max": 2},
        {"key": "dose", "msg": "20. Total amount consumed THIS session (grams/ml)?", "type": "float", "default": 0.3},
        {"key": "drag", "msg": "21. Average puff size? (1: Small/Sips, 2: Normal, 3: Deep lung)", "type": "int", "default": 2, "min": 1, "max": 3},
        {"key": "hold", "msg": "22. Smoke hold time? (1: Instant, 2: 1-2s, 3: 3-5s, 4: >5s)", "type": "int", "default": 2, "min": 1, "max": 4},
        {"key": "duration", "msg": "23. Session duration in minutes?", "type": "float", "default": 10.0},
        {"key": "high_dur", "msg": "24. Compared to friends, your high lasts: (1: Much longer -> 5: Much shorter)", "type": "int", "default": 3, "min": 1, "max": 5},
        {"key": "edible", "msg": "25. Sensitivity to edibles? (1: Too intense, 3: Normal, 5: Barely feel them)", "type": "int", "default": 3, "min": 1, "max": 5},
        {"key": "grogginess", "msg": "26. Next-day grogginess? (1: Always heavy, 3: Sometimes, 5: Never)", "type": "int", "default": 3, "min": 1, "max": 5},
        {"key": "fasting", "msg": "27. Stomach status? (1: Fasted, 2: Normal, 3: High Fat Meal)", "type": "int", "default": 2, "min": 1, "max": 3},
        {"key": "med_list", "custom": "medication_loop"},
        {"key": "post_act", "msg": "29. Activity AFTER consuming? (1: Couch/Sleep, 2: Walking, 3: Heavy Gym)", "type": "int", "default": 1, "min": 1, "max": 3},
    ]

    answers = {}
    i = 0
    while i < len(prompts):
        p = prompts[i]
        if p.get("custom") == "medication_loop":
            answers["med_list"] = run_medication_questionnaire()
            i += 1
            continue

        user_input = input(f"{p['msg']} [Default: {p['default']}]: ").strip()
        if user_input.lower() == 'b':
            if i > 0: i -= 1; print("\n" + "-"*40 + "\n ⏪ Going back...\n" + "-"*40)
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

    timeline, threshold_hour, cleared_hour, rsd, cyp_phenotype, f, c0, base_c0, med_warnings = simulate_thc_clearance(
        answers["gender"], answers["age"], answers["weight"], answers["height"], answers["body_fat"], 
        answers["activity"], answers["weight_trend"], answers["hydration"], answers["diet"],
        answers["years"], answers["days_wk"], answers["grams_day"], answers["break"], answers["source"],
        answers["method"], answers["tobacco"], answers["alcohol"], answers["thc_pct"], answers["cbd"], 
        answers["dose"], answers["drag"], answers["hold"], answers["duration"],
        answers["high_dur"], answers["edible"], answers["grogginess"], answers["fasting"],
        answers["med_list"], answers["post_act"]
    )

    print("\n" + "="*80)
    print(" 🔬 3-COMPARTMENT ELIMINATION ANALYSIS")
    print("="*80)
    
    target_t = threshold_hour if threshold_hour != -1.0 else 48.0
    if target_t <= 6.0: num_steps = 4
    elif target_t <= 24.0: num_steps = 6
    elif target_t <= 72.0: num_steps = 8
    else: num_steps = 10
        
    step_size = target_t / num_steps
    selected_hours = set()
    for k in range(1, num_steps + 1): selected_hours.add(max(0.5, round(k * step_size * 2) / 2.0))
    selected_hours.add(target_t)
    if cleared_hour != -1.0: selected_hours.add(cleared_hour)
        
    for t in sorted(list(selected_hours)):
        data = timeline.get(t)
        if not data: continue
        conc = data["conc"]
        low = data["lower"]
        high = data["upper"]
        status = "🔴 Impaired" if conc >= 3.5 else "🟢 Trace" if conc >= 1.0 else "⚪ Cleared"
        t_str = f"+ {t:.1f} hours" if t < 24 else f"+ {t/24:.1f} days"
        marker = " <== LEGAL THRESHOLD (3.5 ng/ml)" if t == threshold_hour else " <== FULLY CLEARED" if t == cleared_hour else ""
        print(f"{t_str:<15} | {round(conc, 2):<12} | {round(low, 2)} - {round(high, 2):<14} | {status}{marker}")
        if t == threshold_hour and cleared_hour != -1.0 and cleared_hour > threshold_hour: print("-" * 80)
    print("-" * 80)
    
    print("\n📊 Generating Pharmacokinetic Curve with Matplotlib...")
    plot_pharmacokinetics(timeline, threshold_hour, cleared_hour, base_c0)
