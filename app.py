"""
CVRP – Provincia de Cartago  ·  Bloque 03 · Trabajo Grupal · II-1122
Minimización de distancia total con entrega dividida (split delivery)
Métodos: Clarke-Wright heurístico + MIP exacto MTZ (PuLP/CBC)
"""
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import math, time, io

# ──────────────────────────────────────────
#  DATOS BASE
# ──────────────────────────────────────────
CANTONES = {
    0: "CD Cartago",
    1: "Cartago",
    2: "Paraíso",
    3: "La Unión",
    4: "Jiménez",
    5: "Turrialba",
    6: "Alvarado",
    7: "Oreamuno",
    8: "El Guarco",
}

DEMANDA_BASE = {0: 0, 1: 124, 2: 48, 3: 75, 4: 15, 5: 61, 6: 12, 7: 36, 8: 35}

DIST = [
    [0,  0,  9, 10, 34, 34, 20,  6,  6],
    [0,  0,  9, 10, 34, 34, 20,  6,  6],
    [9,  9,  0, 19, 26, 28, 13,  7, 12],
    [10, 10, 19,  0, 45, 43, 29, 14, 11],
    [34, 34, 26, 45,  0, 20, 21, 31, 37],
    [34, 34, 28, 43, 20,  0, 15, 29, 40],
    [20, 20, 13, 29, 21, 15,  0, 14, 25],
    [6,  6,  7, 14, 31, 29, 14,  0, 12],
    [6,  6, 12, 11, 37, 40, 25, 12,  0],
]

# Coordenadas para el mapa (canvas 480×380, y invertido)
COORDS = {
    0: (130, 195),
    1: (148, 215),
    2: (195, 270),
    3: (75,  200),
    4: (270, 245),
    5: (365, 200),
    6: (220, 158),
    7: (175, 118),
    8: (128, 275),
}

COLORS_ROUTES = [
    "#e63946","#2a9d8f","#e9c46a","#f4a261","#6a4c93",
    "#1982c4","#8ac926","#ff595e","#ffca3a","#6a994e",
    "#bc4749","#48cae4","#90be6d","#f9844a","#023e8a",
    "#c77dff","#48a999","#fb8500","#8338ec","#3a86ff",
    "#06d6a0",
]

Q = 24  # capacidad camión

# ──────────────────────────────────────────
#  EXPANSIÓN DE NODOS (split delivery)
# ──────────────────────────────────────────
def build_expanded_nodes(q=Q):
    """Expande cantones con demanda > q en múltiples sub-nodos."""
    nodes = []  # lista de dict {id, canton, load, label}
    for c in range(1, 9):
        dem = DEMANDA_BASE[c]
        n_trips = math.ceil(dem / q)
        remaining = dem
        for t in range(n_trips):
            load = min(q, remaining)
            nodes.append({"id": len(nodes) + 1, "canton": c,
                           "load": load, "label": f"{CANTONES[c]} v{t+1}"})
            remaining -= load
    return nodes

# ──────────────────────────────────────────
#  DISTANCIA ENTRE NODOS EXPANDIDOS
# ──────────────────────────────────────────
def d_exp(ni, nj, exp_nodes):
    """Distancia entre dos nodos expandidos (usan la matriz original por cantón)."""
    ci = exp_nodes[ni - 1]["canton"]
    cj = exp_nodes[nj - 1]["canton"]
    return DIST[ci][cj]

def d_dep(ni, exp_nodes):
    """Distancia del depósito al nodo expandido."""
    ci = exp_nodes[ni - 1]["canton"]
    return DIST[0][ci]

# ──────────────────────────────────────────
#  CLARKE-WRIGHT (sobre nodos expandidos)
# ──────────────────────────────────────────
def clarke_wright(exp_nodes, q=Q):
    n = len(exp_nodes)
    ids = [e["id"] for e in exp_nodes]

    # Ahorros
    savings = []
    for i in ids:
        for j in ids:
            if i < j:
                s = d_dep(i, exp_nodes) + d_dep(j, exp_nodes) - d_exp(i, j, exp_nodes)
                savings.append((s, i, j))
    savings.sort(reverse=True)

    # Rutas iniciales
    routes = {i: [i] for i in ids}
    route_of = {i: i for i in ids}
    load = {i: exp_nodes[i - 1]["load"] for i in ids}

    for s_val, i, j in savings:
        ri = route_of.get(i)
        rj = route_of.get(j)
        if ri is None or rj is None or ri == rj:
            continue
        r_i = routes[ri]
        r_j = routes[rj]
        new_route = None
        if r_i[-1] == i and r_j[0] == j:
            new_route = r_i + r_j
        elif r_i[0] == i and r_j[-1] == j:
            new_route = r_j + r_i
        elif r_i[-1] == i and r_j[-1] == j:
            new_route = r_i + r_j[::-1]
        elif r_i[0] == i and r_j[0] == j:
            new_route = r_i[::-1] + r_j
        else:
            continue
        new_load = load[ri] + load[rj]
        if new_load <= q:
            routes[ri] = new_route
            load[ri] = new_load
            del routes[rj]
            del load[rj]
            for c in new_route:
                route_of[c] = ri

    return list(routes.values()), load

# ──────────────────────────────────────────
#  CÁLCULO DE MÉTRICAS DE RUTA
# ──────────────────────────────────────────
def route_distance_exp(route, exp_nodes):
    if not route:
        return 0
    d = d_dep(route[0], exp_nodes)
    for k in range(len(route) - 1):
        d += d_exp(route[k], route[k + 1], exp_nodes)
    d += d_dep(route[-1], exp_nodes)
    return d

def route_load_exp(route, exp_nodes):
    return sum(exp_nodes[r - 1]["load"] for r in route)

def total_distance_exp(routes, exp_nodes):
    return sum(route_distance_exp(r, exp_nodes) for r in routes)

# ──────────────────────────────────────────
#  MIP EXACTO con PuLP (nodos expandidos)
# ──────────────────────────────────────────
def solve_mip(exp_nodes, q=Q, time_limit=90):
    import pulp
    ids = [e["id"] for e in exp_nodes]
    all_nodes = [0] + ids

    prob = pulp.LpProblem("CVRP_Cartago_Exact", pulp.LpMinimize)

    # Distancias entre todos los nodos
    def d(i, j):
        if i == 0 and j == 0:
            return 0
        if i == 0:
            return d_dep(j, exp_nodes)
        if j == 0:
            return d_dep(i, exp_nodes)
        return d_exp(i, j, exp_nodes)

    # Variables binarias x[i,j]
    x = {(i, j): pulp.LpVariable(f"x_{i}_{j}", cat="Binary")
         for i in all_nodes for j in all_nodes if i != j}
    # MTZ
    u = {i: pulp.LpVariable(f"u_{i}", lowBound=1, upBound=len(ids), cat="Continuous")
         for i in ids}
    # Número de rutas
    K = pulp.LpVariable("K", lowBound=math.ceil(sum(DEMANDA_BASE.values()) / q), cat="Integer")

    # Objetivo
    prob += pulp.lpSum(d(i, j) * x[i, j] for (i, j) in x)

    # Flujo
    for j in ids:
        prob += pulp.lpSum(x[i, j] for i in all_nodes if i != j) == 1
    for i in ids:
        prob += pulp.lpSum(x[i, j] for j in all_nodes if i != j) == 1

    # Rutas desde/hacia depósito
    prob += pulp.lpSum(x[0, j] for j in ids) == K
    prob += pulp.lpSum(x[i, 0] for i in ids) == K

    # MTZ subtour elimination
    M = len(ids)
    for i in ids:
        for j in ids:
            if i != j:
                prob += u[i] - u[j] + M * x[i, j] <= M - 1

    solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=time_limit)
    prob.solve(solver)

    status = pulp.LpStatus[prob.status]
    obj = pulp.value(prob.objective) if prob.objective else None

    routes_mip = []
    if status in ["Optimal", "Feasible"]:
        arcs = [(i, j) for (i, j) in x if pulp.value(x[i, j]) and pulp.value(x[i, j]) > 0.5]
        starts = [j for (i, j) in arcs if i == 0]
        for s in starts:
            route = [s]
            cur = s
            for _ in range(len(ids)):
                nexts = [j for (i, j) in arcs if i == cur and j != 0]
                if not nexts:
                    break
                cur = nexts[0]
                if cur in route:
                    break
                route.append(cur)
            routes_mip.append(route)

    return status, obj, routes_mip

# ──────────────────────────────────────────
#  VISUALIZACIÓN
# ──────────────────────────────────────────
def draw_map(routes, exp_nodes, title="CVRP · Cartago"):
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set_facecolor("#0d1b2a")
    fig.patch.set_facecolor("#0d1b2a")
    ax.set_xlim(30, 450)
    ax.set_ylim(70, 340)
    ax.invert_yaxis()
    ax.axis("off")

    canton_coords = COORDS  # índice por cantón base

    # Dibujar rutas
    for idx, route in enumerate(routes):
        color = COLORS_ROUTES[idx % len(COLORS_ROUTES)]
        cantons_route = [0] + [exp_nodes[r - 1]["canton"] for r in route] + [0]
        xs = [canton_coords[c][0] for c in cantons_route]
        ys = [canton_coords[c][1] for c in cantons_route]

        # Línea con offset ligero para rutas al mismo cantón
        ax.plot(xs, ys, "-", color=color, linewidth=1.8, alpha=0.7, zorder=2)
        # Flechas
        for k in range(len(xs) - 1):
            dx = xs[k + 1] - xs[k]
            dy = ys[k + 1] - ys[k]
            ax.annotate("", xy=(xs[k + 1], ys[k + 1]), xytext=(xs[k], ys[k]),
                        arrowprops=dict(arrowstyle="-|>", color=color,
                                        lw=1.4, mutation_scale=10),
                        zorder=3)

    # Nodos
    for nid, (xc, yc) in canton_coords.items():
        is_depot = nid == 0
        color = "#ffd700" if is_depot else "#e0e0e0"
        size = 180 if is_depot else 100
        ax.scatter(xc, yc, s=size, color=color, zorder=5,
                   edgecolors="#0d1b2a", linewidths=1.5)
        label = CANTONES[nid]
        offset_x, offset_y = (6, -14) if nid != 0 else (7, -16)
        txt = ax.text(xc + offset_x, yc + offset_y, label,
                      fontsize=7.5, color="white", zorder=6,
                      fontweight="bold" if is_depot else "normal")
        txt.set_path_effects([pe.withStroke(linewidth=2, foreground="#0d1b2a")])

    # Leyenda de rutas
    handles = []
    for idx, route in enumerate(routes):
        ld = route_load_exp(route, exp_nodes)
        km = route_distance_exp(route, exp_nodes)
        label = f"R{idx+1}: {ld}pal · {km}km"
        handles.append(plt.Line2D([0], [0], color=COLORS_ROUTES[idx % len(COLORS_ROUTES)],
                                   lw=2, label=label))
    ax.legend(handles=handles, loc="lower right", fontsize=6.5,
              facecolor="#0d1b2a", labelcolor="white",
              edgecolor="#333", framealpha=0.95, ncol=2)

    total_km = total_distance_exp(routes, exp_nodes)
    ax.set_title(f"{title}\nDistancia total: {total_km:.0f} km · {len(routes)} viajes",
                 color="white", fontsize=11, fontweight="bold", pad=8)
    plt.tight_layout()
    return fig

# ──────────────────────────────────────────
#  STREAMLIT APP
# ──────────────────────────────────────────
st.set_page_config(
    page_title="CVRP · Cartago",
    page_icon="🚛",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Syne:wght@400;700;800&family=DM+Sans:wght@300;400;500&display=swap');
html,body,[class*="css"]{font-family:'DM Sans',sans-serif;}
h1,h2,h3{font-family:'Syne',sans-serif!important;}
.kpi{background:linear-gradient(135deg,#16213e,#0f3460);border:1px solid #1a4a8a;
     border-radius:12px;padding:14px 18px;margin-bottom:8px;}
.kpi-val{font-size:1.9rem;font-weight:800;color:#e2d9f3;line-height:1;}
.kpi-lbl{font-size:0.72rem;color:#90a0b0;text-transform:uppercase;letter-spacing:.08em;margin-top:4px;}
.stButton>button{background:linear-gradient(90deg,#e63946,#c1121f);color:white;
                 border:none;border-radius:8px;font-family:'Syne',sans-serif;
                 font-weight:700;letter-spacing:.05em;}
</style>
""", unsafe_allow_html=True)

# ── HEADER ──
st.markdown("""
<div style="background:linear-gradient(135deg,#0d1b2a,#1b263b);border-radius:16px;
            padding:26px 30px;margin-bottom:22px;border-left:5px solid #e63946;">
  <div style="font-size:.7rem;color:#f4a261;font-weight:700;
              letter-spacing:.15em;text-transform:uppercase;margin-bottom:6px;">
    Bloque 03 · Trabajo Grupal · II-1122
  </div>
  <h1 style="margin:0;color:white;font-size:2rem;">🚛 CVRP · Provincia de Cartago</h1>
  <p style="color:#a0aec0;margin:7px 0 0;font-size:.88rem;">
    Minimización de distancia total · Entrega dividida (split delivery) · 
    8 cantones · 406 pallets/semana · Q = 24 pallets/camión
  </p>
</div>
""", unsafe_allow_html=True)

# ── SIDEBAR ──
with st.sidebar:
    st.markdown("### ⚙️ Parámetros")
    q_cap = st.slider("Capacidad camión (pallets)", 12, 48, 24, 4)
    metodo = st.radio("Método", ["Clarke-Wright (Heurística)", "MIP Exacto (PuLP/CBC)"])

    if "MIP" in metodo:
        tl = st.slider("Tiempo límite solver (s)", 30, 300, 90, 30)
    else:
        tl = 90

    st.markdown("---")
    st.markdown("### 📊 Demanda por cantón")
    df_dem = pd.DataFrame([{"Nodo": k, "Cantón": CANTONES[k],
                             "Demanda (pal)": DEMANDA_BASE[k],
                             "Viajes min": math.ceil(DEMANDA_BASE[k] / q_cap) if k > 0 else 0}
                            for k in range(9) if k > 0])
    st.dataframe(df_dem, hide_index=True, use_container_width=True)
    total_dem = sum(v for k, v in DEMANDA_BASE.items() if k > 0)
    total_viajes = sum(math.ceil(v / q_cap) for k, v in DEMANDA_BASE.items() if k > 0)
    st.info(f"**Demanda total:** {total_dem} pal/sem  \n"
            f"**Viajes mínimos:** {total_viajes}  \n"
            f"**Flota mínima:** ⌈{total_dem}/{q_cap}⌉ = **{math.ceil(total_dem/q_cap)}**")

# ── TABS ──
tab1, tab2, tab3, tab4 = st.tabs(["🗺️ Solución & Mapa", "📋 Detalle de Viajes",
                                   "📐 Datos del Problema", "📖 Formulación MIP"])

# ────────── TAB 1 ──────────
with tab1:
    c_btn, _ = st.columns([1, 3])
    with c_btn:
        run = st.button("▶ Resolver CVRP", use_container_width=True)

    if run:
        exp_nodes = build_expanded_nodes(q_cap)
        with st.spinner(f"Resolviendo con {metodo}..."):
            t0 = time.time()
            if "Clarke-Wright" in metodo:
                routes_r, loads_r = clarke_wright(exp_nodes, q_cap)
                status_lbl = "Clarke-Wright (Heurística)"
                obj_val = total_distance_exp(routes_r, exp_nodes)
            else:
                st.info(f"Ejecutando MIP exacto (límite: {tl}s)...")
                status, obj, routes_mip = solve_mip(exp_nodes, q_cap, tl)
                if routes_mip:
                    routes_r = routes_mip
                    status_lbl = f"MIP Exacto · {status}"
                    obj_val = obj or total_distance_exp(routes_r, exp_nodes)
                else:
                    routes_r, loads_r = clarke_wright(exp_nodes, q_cap)
                    status_lbl = f"MIP sin sol. factible (usó CW) · {status}"
                    obj_val = total_distance_exp(routes_r, exp_nodes)
            elapsed = time.time() - t0

        st.session_state.update({
            "routes": routes_r, "exp_nodes": exp_nodes,
            "obj": obj_val, "status": status_lbl, "elapsed": elapsed
        })

    if "routes" in st.session_state:
        routes_r = st.session_state["routes"]
        exp_nodes = st.session_state["exp_nodes"]
        obj_val = st.session_state["obj"]
        status_lbl = st.session_state["status"]
        elapsed = st.session_state["elapsed"]

        n_viajes = len(routes_r)
        total_km = total_distance_exp(routes_r, exp_nodes)
        total_dem_v = sum(e["load"] for e in exp_nodes)

        c1, c2, c3, c4 = st.columns(4)
        for col, val, lbl in [
            (c1, f"{total_km:.0f} km", "Distancia Total"),
            (c2, str(n_viajes), "Viajes / Rutas"),
            (c3, str(math.ceil(sum(DEMANDA_BASE.values()) / q_cap)), "Flota mínima"),
            (c4, f"{elapsed:.1f}s", "Tiempo CPU"),
        ]:
            col.markdown(f'<div class="kpi"><div class="kpi-val">{val}</div>'
                         f'<div class="kpi-lbl">{lbl}</div></div>', unsafe_allow_html=True)

        st.caption(f"**Método:** {status_lbl}")

        fig = draw_map(routes_r, exp_nodes)
        st.pyplot(fig, use_container_width=True)

        # Botón descarga PNG
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        st.download_button("⬇ Descargar mapa PNG", buf, "cvrp_cartago_rutas.png", "image/png")
    else:
        st.info("Presioná **▶ Resolver CVRP** para generar las rutas óptimas.")

# ────────── TAB 2 ──────────
with tab2:
    st.markdown("### Detalle de viajes")
    if "routes" in st.session_state:
        routes_r = st.session_state["routes"]
        exp_nodes = st.session_state["exp_nodes"]
        q_used = q_cap

        rows = []
        for idx, route in enumerate(routes_r):
            stops_raw = [exp_nodes[r - 1]["label"] for r in route]
            # Simplificar: mostrar cantones únicos visitados
            cantons_visited = list(dict.fromkeys(exp_nodes[r - 1]["canton"] for r in route))
            stops_str = " → ".join(["CD Cartago"] +
                                    [CANTONES[c] for c in cantons_visited] + ["CD Cartago"])
            ld = route_load_exp(route, exp_nodes)
            km = route_distance_exp(route, exp_nodes)
            rows.append({
                "Viaje": f"V{idx + 1}",
                "Ruta": stops_str,
                "Carga (pal)": ld,
                "Dist. (km)": round(km, 1),
                "% capacidad": f"{ld / q_used * 100:.0f}%",
            })

        df_r = pd.DataFrame(rows)
        st.dataframe(df_r, hide_index=True, use_container_width=True)

        avg_km = df_r["Dist. (km)"].mean()
        avg_ld = df_r["Carga (pal)"].mean()
        total_km2 = df_r["Dist. (km)"].sum()
        st.markdown(f"""
        **Total distancia:** {total_km2:.1f} km &nbsp;|&nbsp;
        **Promedio km/viaje:** {avg_km:.1f} km &nbsp;|&nbsp;
        **Carga promedio:** {avg_ld:.1f} pal &nbsp;|&nbsp;
        **Eficiencia promedio:** {avg_ld / q_used * 100:.0f}%
        """)

        # Verificación
        covered_cantons = set(exp_nodes[r - 1]["canton"]
                               for route in routes_r for r in route)
        all_cantons = set(range(1, 9))
        if covered_cantons == all_cantons:
            st.success("✅ Los 8 cantones están cubiertos en la solución.")
        else:
            missing = all_cantons - covered_cantons
            st.warning(f"⚠️ Cantones faltantes: {[CANTONES[c] for c in missing]}")
    else:
        st.info("Ejecutá la solución primero.")

# ────────── TAB 3 ──────────
with tab3:
    st.markdown("### Demanda por cantón")
    df_d = pd.DataFrame([
        {"Nodo": k, "Cantón": CANTONES[k], "Imperial": [62,24,37,7,31,6,18,17][k-1],
         "Pilsen": [31,12,19,4,15,3,9,9][k-1], "Tropical": [31,12,19,4,15,3,9,9][k-1],
         "Demanda total": DEMANDA_BASE[k]}
        for k in range(1, 9)
    ])
    st.dataframe(df_d, hide_index=True, use_container_width=True)

    st.markdown("### Matriz de distancias por carretera (km)")
    df_dist = pd.DataFrame(DIST,
                           index=[f"{i}·{CANTONES[i]}" for i in range(9)],
                           columns=[str(i) for i in range(9)])
    st.dataframe(df_dist.style.background_gradient(cmap="OrRd", axis=None)
                              .format("{:.0f}"), use_container_width=True)

    st.markdown("### Tabla de ahorros Clarke-Wright")
    sav_rows = []
    for i in range(1, 9):
        for j in range(1, 9):
            if i < j:
                s = DIST[0][i] + DIST[0][j] - DIST[i][j]
                sav_rows.append({"i": CANTONES[i], "j": CANTONES[j],
                                  "d(0,i)": DIST[0][i], "d(0,j)": DIST[0][j],
                                  "d(i,j)": DIST[i][j], "Ahorro s(i,j)": s})
    df_sav = pd.DataFrame(sav_rows).sort_values("Ahorro s(i,j)", ascending=False)
    st.dataframe(df_sav, hide_index=True, use_container_width=True)

# ────────── TAB 4 ──────────
with tab4:
    st.markdown("### Formulación MIP del CVRP (Miller-Tucker-Zemlin)")
    st.markdown(r"""
**Conjuntos y parámetros:**
- $N = \{0,1,\ldots,n\}$ : nodos expandidos (0 = depósito CD Cartago; $n$ = total sub-nodos)
- $C \subset N$ : sub-nodos clientes tras expansión por split delivery
- $d_{ij}$ : distancia por carretera (km) según la matriz original
- $q_i$ : carga asignada al sub-nodo $i$ ($\leq Q$)  
- $Q = 24$ : capacidad de cada camión

**Variables:**

$$x_{ij} \in \{0,1\} \quad \text{(arco $(i,j)$ recorrido)}$$
$$u_i \geq 0 \quad \text{(posición en ruta, anti-subtour MTZ)}$$
$$K \in \mathbb{Z}^+ \quad \text{(número de rutas)}$$

**Función objetivo:**

$$\min \quad Z = \sum_{i \in N}\sum_{\substack{j \in N \\ j \neq i}} d_{ij}\, x_{ij}$$

**Restricciones:**

$$\sum_{i \in N,\, i \neq j} x_{ij} = 1 \qquad \forall j \in C \tag{1 · Cobertura entrada}$$

$$\sum_{j \in N,\, j \neq i} x_{ij} = 1 \qquad \forall i \in C \tag{2 · Cobertura salida}$$

$$\sum_{j \in C} x_{0j} = K \tag{3 · Salidas del depósito}$$

$$\sum_{i \in C} x_{i0} = K \tag{4 · Regresos al depósito}$$

$$u_i - u_j + |C|\cdot x_{ij} \leq |C|-1 \qquad \forall i,j \in C,\; i \neq j \tag{5 · MTZ anti-subtour}$$

$$K \geq \left\lceil \frac{\displaystyle\sum_{i \in C} q_i}{Q} \right\rceil = \left\lceil \frac{406}{24} \right\rceil = 17 \tag{6 · Flota mínima}$$

$$x_{ij} \in \{0,1\}, \quad u_i \geq 0, \quad K \in \mathbb{Z}^+ \tag{7 · Dominio}$$

---
**Nota:** Con entrega dividida se generan **21 sub-nodos**. El MIP MTZ sobre 21 nodos 
puede resolverse con CBC en ≈ 30–90 s a optimalidad. Para la heurística Clarke-Wright,
se aplica el mismo algoritmo de ahorros sobre los sub-nodos expandidos.
    """)

# ── FOOTER ──
st.markdown("""
<hr style="border-color:#2d3748;margin-top:36px;">
<p style="text-align:center;color:#718096;font-size:.76rem;">
Prof. David Benavides · UCR Sede Alajuela · I-2026 · Bloque 03 – Trabajo Grupal · II-1122
</p>
""", unsafe_allow_html=True)
