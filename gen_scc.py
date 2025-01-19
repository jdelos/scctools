import numpy as np
import sympy as sp

class GenericSwitchedCapacitor:
    def __init__(self, ArchDesc, **kwargs):
        # Validate and extract architecture matrices
        self.inc_caps = ArchDesc.get('Acaps')
        self.inc_switches = ArchDesc.get('Asw')
        self.sw_activation_matrix = ArchDesc.get('Asw_act')

        if self.inc_caps is None or self.inc_switches is None or self.sw_activation_matrix is None:
            raise ValueError("Missing required fields in ArchDesc: 'Acaps', 'Asw', or 'Asw_act'")

        # Basic attributes
        self.n_phases = self.sw_activation_matrix.shape[0]
        self.n_nodes = self.inc_caps.shape[0]
        self.n_caps = self.inc_caps.shape[1]
        self.n_switches = self.inc_switches.shape[1]
        self.duty = kwargs.get('Duty', sp.symbols(f'D1:{self.n_phases}'))
        self.mode = kwargs.get('Mode', 0)
        self.in_node = kwargs.get('InNode', 1)
        self.fsw_op = 1  # Normalized switching frequency

        # Initialize symbolic parameters
        self.ron_switches = sp.symbols(f'Ron1:{self.n_switches + 1}')
        self.caps = sp.symbols(f'C1:{self.n_caps + 1}')
        self.esr_caps = sp.symbols(f'Resr1:{self.n_caps + 1}')
        self.phase = [{} for _ in range(self.n_phases)]

        # Generate activation matrices
        self.inc_sw_act = [self.inc_switches[:, self.sw_activation_matrix[j, :].astype(bool)] for j in range(self.n_phases)]
        self.incidence_matrix = np.concatenate([self.inc_caps] + self.inc_sw_act, axis=1)

        # Generate initial calculations
        self.gen_k()
        self.caps_voltage_ratio()
        self.switch_voltage_ratio()

    def gen_k(self):
        # Compute slow and fast switching limit resistances (simplified)
        self.k_ssl = 1 / (2 * self.fsw_op) * sum([sp.Matrix(self.inc_caps).norm() for _ in range(self.n_phases)])
        self.k_fsl = sum([sp.Matrix(self.inc_switches).norm() for _ in range(self.n_phases)])
        self.k_factors = (self.k_ssl**2 + self.k_fsl**2)**0.5

    def caps_voltage_ratio(self):
        # Simplified placeholder for voltage normalization at capacitors
        B = np.vstack([self.inc_caps for _ in range(self.n_phases)])
        self.v_caps_norm = -np.linalg.pinv(B[:, 1:]) @ B[:, 0]

    def switch_voltage_ratio(self):
        # Simplified placeholder for voltage normalization at switches
        sw_volt = np.zeros((self.n_phases, self.n_switches))
        self.v_sw_norm = np.max(np.abs(sw_volt), axis=0)

    def Rssl(self, Fsw, Cx):
        # Compute slow switching limit resistance
        if isinstance(Cx, (int, float)):
            Cx = [Cx] * self.n_caps
        return [float(1 / (2 * Fsw * c)) for c in Cx]

    def Rfsl(self, Ron):
        # Compute fast switching limit resistance
        if isinstance(Ron, (int, float)):
            Ron = [Ron] * self.n_switches
        return [float(r) for r in Ron]

    def Resr(self, Cesr):
        # Compute ESR-related resistance
        if isinstance(Cesr, (int, float)):
            Cesr = [Cesr] * self.n_caps
        return [float(r) for r in Cesr]

    def Rscc(self, Fsw, Cx, Ron, Cesr):
        # Combined resistance computation
        rssl = self.Rssl(Fsw, Cx)
        rfsl = self.Rfsl(Ron)
        resr = self.Resr(Cesr)
        return [((rssl[i] ** 2 + (rfsl[i] + resr[i]) ** 2) ** 0.5) for i in range(len(rssl))]

# Example usage
ArchDesc = {
    'Acaps': np.array([[1, -1], [-1, 1]]),
    'Asw': np.array([[1, 0], [0, 1]]),
    'Asw_act': np.array([[1, 0], [0, 1]])
}

converter = GenericSwitchedCapacitor(ArchDesc, Mode=0, InNode=1, Duty=[0.5, 0.5])
print("SSL Resistance:", converter.Rssl(100e3, [1e-6, 1e-6]))
print("FSL Resistance:", converter.Rfsl([0.01, 0.01]))
print("SCC Resistance:", converter.Rscc(100e3, [1e-6, 1e-6], [0.01, 0.01], [0.001, 0.001]))
