# 🚛 CVRP – Provincia de Cartago

**Bloque 03 · Trabajo Grupal · Curso II-1122 – Programación Entera Mixta con Software**  
Prof. David Benavides · UCR Sede Alajuela · I-2026

---

## Descripción del problema

Resolver el **Capacitated Vehicle Routing Problem (CVRP)** para la distribución de productos
cerveceros (Imperial / Pilsen / Tropical) desde el Centro de Distribución (CD) en Cartago hacia
los 8 cantones de la provincia.

| Parámetro | Valor |
|---|---|
| Nodos | 9 (1 depósito + 8 cantones) |
| Demanda total | 406 pallets/semana |
| Capacidad camión | 24 pallets |
| Flota mínima | ⌈406/24⌉ = 17 camiones |
| Objetivo | Minimizar distancia total (km) |

## Estructura del repositorio

```
cvrp_cartago/
├── app.py               # Aplicación Streamlit principal
├── requirements.txt     # Dependencias Python
├── README.md            # Este archivo
└── .streamlit/
    └── config.toml      # Configuración de tema Streamlit
```

## Métodos implementados

### 1. Heurística Clarke-Wright (Savings Algorithm)
- Calcula ahorros `s(i,j) = d(0,i) + d(0,j) - d(i,j)` para todos los pares de clientes
- Fusiona rutas en orden descendente de ahorro respetando la capacidad Q = 24
- Complejidad: O(n² log n)

### 2. MIP Exacto con formulación MTZ (PuLP + CBC)
**Función objetivo:**
```
min  Σ_{i,j} d_ij · x_ij
```
**Restricciones principales:**
- Flujo de entrada/salida por cantón (= 1)
- Eliminación de sub-tours (Miller-Tucker-Zemlin)
- Capacidad de carga por ruta
- Flota mínima

## Instalación local

```bash
git clone https://github.com/<tu-usuario>/cvrp_cartago.git
cd cvrp_cartago
pip install -r requirements.txt
streamlit run app.py
```

## Despliegue en Streamlit Cloud

1. Sube el repositorio a GitHub (público o privado)
2. Ve a [share.streamlit.io](https://share.streamlit.io)
3. Conecta el repositorio y selecciona `app.py` como archivo principal
4. Click en **Deploy** — listo en ~2 minutos

## Datos

### Demanda por cantón (pallets/semana)

| Nodo | Cantón | Imperial | Pilsen | Tropical | Total |
|---|---|---|---|---|---|
| 0 | CD Cartago | - | - | - | — |
| 1 | Cartago | 62 | 31 | 31 | **124** |
| 2 | Paraíso | 24 | 12 | 12 | **48** |
| 3 | La Unión | 37 | 19 | 19 | **75** |
| 4 | Jiménez | 7 | 4 | 4 | **15** |
| 5 | Turrialba | 31 | 15 | 15 | **61** |
| 6 | Alvarado | 6 | 3 | 3 | **12** |
| 7 | Oreamuno | 18 | 9 | 9 | **36** |
| 8 | El Guarco | 17 | 9 | 9 | **35** |

### Matriz de distancias por carretera (km)

|  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| **0** | 0 | 0 | 9 | 10 | 34 | 34 | 20 | 6 | 6 |
| **1** | 0 | 0 | 9 | 10 | 34 | 34 | 20 | 6 | 6 |
| **2** | 9 | 9 | 0 | 19 | 26 | 28 | 13 | 7 | 12 |
| **3** | 10 | 10 | 19 | 0 | 45 | 43 | 29 | 14 | 11 |
| **4** | 34 | 34 | 26 | 45 | 0 | 20 | 21 | 31 | 37 |
| **5** | 34 | 34 | 28 | 43 | 20 | 0 | 15 | 29 | 40 |
| **6** | 20 | 20 | 13 | 29 | 21 | 15 | 0 | 14 | 25 |
| **7** | 6 | 6 | 7 | 14 | 31 | 29 | 14 | 0 | 12 |
| **8** | 6 | 6 | 12 | 11 | 37 | 40 | 25 | 12 | 0 |

## Integrantes del grupo

- _(Agregar nombres aquí)_

---
*Curso II-1122 · UCR Sede Alajuela · I-2026*
