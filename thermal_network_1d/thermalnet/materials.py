"""材料物性ライブラリ(代表値)。

k   : 熱伝導率 [W/(m K)]
rho : 密度 [kg/m^3]
cp  : 比熱 [J/(kg K)]

注意: ここの値はあくまで「教材用の代表値」。実設計ではメーカーの
データシート値を使うこと。特に樹脂・TIM・基板は品種で 2-3 倍変わる。
"""

MATERIALS = {
    # 金属
    "aluminum":       {"k": 167.0, "rho": 2700.0, "cp": 900.0},
    "aluminum_6063":  {"k": 200.0, "rho": 2700.0, "cp": 900.0},
    "copper":         {"k": 390.0, "rho": 8933.0, "cp": 385.0},
    "steel":          {"k": 50.0,  "rho": 7850.0, "cp": 470.0},
    "stainless_304":  {"k": 15.0,  "rho": 7900.0, "cp": 477.0},
    "solder_sac":     {"k": 58.0,  "rho": 7400.0, "cp": 230.0},
    # 半導体・セラミック
    "silicon":        {"k": 120.0, "rho": 2330.0, "cp": 700.0},
    "sic":            {"k": 370.0, "rho": 3210.0, "cp": 690.0},
    "alumina":        {"k": 30.0,  "rho": 3900.0, "cp": 880.0},
    "aln":            {"k": 170.0, "rho": 3260.0, "cp": 740.0},
    # 基板・樹脂
    "fr4":            {"k": 0.35,  "rho": 1900.0, "cp": 1150.0},
    "abs":            {"k": 0.20,  "rho": 1050.0, "cp": 1400.0},
    "epoxy_mold":     {"k": 0.80,  "rho": 1900.0, "cp": 900.0},
    # 建築
    "gypsum_board":   {"k": 0.22,  "rho": 800.0,  "cp": 1000.0},
    "glass_wool":     {"k": 0.040, "rho": 20.0,   "cp": 840.0},
    "plywood":        {"k": 0.15,  "rho": 550.0,  "cp": 1600.0},
    "glass":          {"k": 1.0,   "rho": 2500.0, "cp": 750.0},
    "concrete":       {"k": 1.6,   "rho": 2300.0, "cp": 880.0},
    "wood_stud":      {"k": 0.13,  "rho": 500.0,  "cp": 1600.0},
    # 流体(静止層としての k、および熱容量流量の cp に使う)
    "air":            {"k": 0.026, "rho": 1.16,   "cp": 1007.0},
    "water":          {"k": 0.60,  "rho": 997.0,  "cp": 4180.0},
    "glycol_50":      {"k": 0.38,  "rho": 1060.0, "cp": 3300.0},
    "oil":            {"k": 0.14,  "rho": 860.0,  "cp": 2000.0},
    # インタフェース材(k のみ。実際は厚みと接触抵抗で決まる)
    "thermal_grease": {"k": 3.0,   "rho": 2500.0, "cp": 1000.0},
    "thermal_pad":    {"k": 3.0,   "rho": 2000.0, "cp": 1000.0},
}


def get(name: str) -> dict:
    """材料名から物性を引く。未登録ならエラーで候補を出す。"""
    key = str(name).strip().lower()
    if key not in MATERIALS:
        raise KeyError(
            f"材料 '{name}' は未登録です。使えるのは: "
            f"{', '.join(sorted(MATERIALS))}\n"
            "  モデルファイルの materials: セクションで自分で定義もできます。")
    return MATERIALS[key]


def property_of(name: str, prop: str) -> float:
    mat = get(name)
    if prop not in mat:
        raise KeyError(f"材料 '{name}' に物性 '{prop}' がありません")
    return mat[prop]
