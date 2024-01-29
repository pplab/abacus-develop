import numpy as np
import scipy.linalg as la
import abacus2qo.components.safe_guard as sg

class toQO_Calculator:
    """python-end the Quasiatomic orbital (QO) analysis
    """
    sg_: sg.toQO_SafeGuard

    def __init__(self) -> None:
        self.sg_ = sg.toQO_SafeGuard()

    def folding_HR(self, matrices_R: list, supercells: list, kpoint: np.ndarray) -> np.ndarray:
        """calculate fold matrix from R to k-space

        Args:
            matrices_R (list): matrices in R-space
            supercells (list): corresponding supercells direct vectors n1n2n3
            kpoint (np.ndarray): one specific kpoint direct coordinates

        Returns:
            np.ndarray: matrix in k-space
        """
        matrix_k = np.zeros_like(matrices_R[0])
        for iR in range(len(supercells)):
            arg = np.exp(1j * supercells[iR] @ kpoint * 2 * np.pi)
            matrix_k += arg * matrices_R[iR]
        return matrix_k
    
    def unfolding_Hk(self, matrices_k: list, kpoints: list, supercell: np.ndarray) -> np.ndarray:
        """calculate unfold matrix from k to R-space

        Args:
            matrices_k (list): matrices in k-space
            kpoints (list): corresponding equivalent kpoints, arranged in the same order as matrices_k
            supercell (np.ndarray): one specific supercell direct coordinates

        Returns:
            np.ndarray: matrix in R-space
        """

        matrix_R = np.zeros_like(matrices_k[0])
        for ik in range(len(kpoints)):
            kpoint = kpoints[ik]
            arg = np.exp(-1j * kpoint @ supercell * 2 * np.pi)
            matrix_R += arg * matrices_k[ik]/np.sqrt(len(kpoints))
        return matrix_R
    
    def projto_nao(self, sk: np.ndarray, saok: np.ndarray) -> np.ndarray:
        """calculate anything in nao representation with the nao projector

        Args:
            sk (np.ndarray): nao overlap matrix in k-space
            saok (np.ndarray): ao-nao overlap matrix in k-space

        Returns:
            np.ndarray: anything in nao representation
        """
        print("Calculate Atomic Orbital (AO) in NAO representation.")
        psi_chi = np.linalg.solve(sk, saok.conj().T)
        return psi_chi

    def projto_eigstate(self, psi_lcao: np.ndarray, saok: np.ndarray) -> np.ndarray:
        """calculate any states projected onto eigenstates of Hamiltonian in k-space represented by NAO

        Args:
            psi_lcao (np.ndarray): states in LCAO representation
            saok (np.ndarray): ao-nao overlap matrix in k-space

        Returns:
            np.ndarray: states projected onto eigenstates of Hamiltonian in k-space represented by NAO
        """
        psi_chi_para = psi_lcao @ psi_lcao.conj().T @ saok.conj().T
        return psi_chi_para
    
    def canonical_orthogonalization(self, psi_chi_orth: np.ndarray, sk: np.ndarray, m: int) -> np.ndarray:
        """this is a general canonical orthogonalization procedure for any subspace

        Args:
            psi_chi_orth (np.ndarray): states in one certain representation, here is nao
            sk (np.ndarray): basis function overlap matrix
            m (int): number of states retained

        Returns:
            np.ndarray: states orthogonalized
        """
        print("Perform canonical orthogonalization for newly appended subspace from AO introduction.")
        wk = psi_chi_orth.conj().T @ sk @ psi_chi_orth
        yk_k, vk_k = la.eigh(wk)
        print("Eigenvalues of W(k) are: \n", yk_k)
        # truncate m largest eigenvalues and eigenvectors
        yk_k = yk_k[-m:]
        vk_k = vk_k[:, -m:]
        print("Get ", m, " largest eigenvalues and eigenvectors. Selected eigenvalues are: \n", yk_k)
        # calculate psi_complem
        yk_k = np.diag([np.sqrt(1/yk_k[i]) for i in range(m)])
        psi_complem = psi_chi_orth @ vk_k @ yk_k
        print("Check orthogonalities between states in extended space.")
        self.sg_.check_eigenstates(psi_complem, sk)
        return psi_complem
    
    def merge_space(self, psi_lcao: np.ndarray, psi_complem: np.ndarray, hk: np.ndarray, sk: np.ndarray) -> np.ndarray:
        """calculate psi_exten

        psi_exten = [psi_lcao, psi_complem]

        Args:
            psi_lcao (np.ndarray): states in LCAO representation
            psi_complem (np.ndarray): states in complementary representation
            hk (np.ndarray): Hamiltonian in k-space
            sk (np.ndarray): overlap in k-space

        Returns:
            np.ndarray: extended states in k-space
        """
        print("Combine occupied states and constructed virtual states to get extended wavefunction"
             +"\n(empty states are appended)."
             +"\nShape of occupied states is: ", psi_lcao.shape, "\nShape of empty states is: ", psi_complem.shape)
        psi_exten = np.concatenate((psi_lcao, psi_complem), axis=1)
        print("Shape of extended states is: ", psi_exten.shape)
        # check orthogonality
        print("Check orthogonality of psi_exten.")
        sk_exten = psi_exten.conj().T @ sk @ psi_exten
        # if sk_exten is not identity matrix?
        error = self.sg_.check_identity(sk_exten)
        while error > 1e-6:
            print("Error of sk_exten not being identity is: ", error, ", unacceptable.")
            exit()
            
        # if hk_exten is not diagonal in supspace psi_lcao?
        hk_exten = psi_exten.conj().T @ hk @ psi_exten
        self.sg_.check_diagonal(hk_exten[:psi_lcao.shape[1], :psi_lcao.shape[1]])
        # get extended energy spectrum
        print("Get extended energy spectrum.")
        eigvals_exten, eigvecs_exten = la.eigh(hk_exten, sk_exten)
        print("Eigenvalues of H_exten(k) are: \n", eigvals_exten)
        eigvals, eigvecs = la.eigh(hk, sk)
        print("Comparatively eigenvalues of H(k) are: \n", eigvals)
        return psi_exten

    def calculate_qo(self, saok: np.ndarray, psi_exten: np.ndarray, sk: np.ndarray) -> np.ndarray:
        """calculate qo represented by NAO in kspace, return with nrows = nchi and ncols = nphi

        Args:
            saok (np.ndarray): ao-nao overlap matrix in k-space
            psi_exten (np.ndarray): extended states in k-space
            sk (np.ndarray): overlap in k-space

        Returns:    
            np.ndarray: qo represented by NAO in kspace
        """
        print("Calculate QO.")
        qo = psi_exten @ psi_exten.conj().T @ saok.conj().T
             # this saok is overlap between AO and NAO, line is AO, column is NAO
        # then normalize qo
        for i in range(qo.shape[1]):
            qo[:, i] = qo[:, i] / np.sqrt(qo[:, i].conj().T @ sk @ qo[:, i])
            print("QO Normalization: after, norm of the ", i, "-th QO is: %10.8f"%((qo[:, i].conj().T @ sk @ qo[:, i]).real))
        return qo

    def calculate_hqok(self, qo: np.ndarray, hk: np.ndarray, sk: np.ndarray) -> np.ndarray:
        """calculate hqok

        Args:
            qo (np.ndarray): QO in k-space
            hk (np.ndarray): Hamiltonian represented by QO in k-space
            sk (np.ndarray): QO overlap in k-space

        Returns:
            np.ndarray: hqok
        """
        print("Calculate hamiltonian matrix in QO basis in k-space.")
        hqok = qo.conj().T @ hk @ qo
        sqok = qo.conj().T @ sk @ qo # this is overlap matrix in QO basis

        skip_diagonalization = False
        eigval_s, eigvec_s = la.eigh(sqok)
        print("Eigenvalues of overlap of in basis Sqo(k) are: \n", eigval_s)
        if abs(eigval_s[0]) < 1e-6:
            print("WARNING: BAD CONDITION NUMBER OF OVERLAP MATRIX IN QO BASIS, THIS WILL CAUSE PROBLEMS IN DIAGONALIZATION.")
            return np.zeros_like(hqok), np.zeros_like(sqok)
        
        else:
            eigvals_qo, eigvecs_qo = la.eigh(hqok, sqok)
            print("Eigenvalues of Hamiltonian in QO basis Hqo(k) are: \n", eigvals_qo)
            
            eigvals_nao, eigvecs_nao = la.eigh(hk, sk)
            print("Eigenvalues of Hamiltonian in NAO basis H(k) are: \n", eigvals_nao)

            return hqok, sqok
