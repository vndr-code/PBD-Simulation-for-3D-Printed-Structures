# spring.py
import numpy as np
import config

class Spring:
    def __init__(self, particle1, particle2, rest_length, stiffness, plasticity, yield_ratio, secondary=False):
        """
        Represents a spring connecting two particles by their indices.

        Parameters:
            particle1 (int): Index of the first particle.
            particle2 (int): Index of the second particle.
            rest_length (float): The rest length of the spring.
            stiffness (float): The spring stiffness constant.
            plasticity (float): The plasticity constant (α).
            yield_ratio (float): The yield ratio (γ).
        """
        self.particle1 = particle1
        self.particle2 = particle2
        self.rest_length = rest_length
        self.stiffness = stiffness
        self.plasticity = plasticity
        self.yield_ratio = yield_ratio
        self.secondary = secondary

    
    def project(self, particles):
        """
        Projects the spring constraint by adjusting the positions of the two connected particles.
        For secondary springs, the effective stiffness, plasticity, and yield_ratio are modified based
        on the (pre-assigned) age of the particles. The idea is that older (bottom) layers become stiffer
        and less plastic.
        """
        p1 = particles[self.particle1]
        p2 = particles[self.particle2]
        diff = p2.pos - p1.pos
        current_length = np.linalg.norm(diff)
        if current_length < 1e-6:
            return
        
        # Set effective parameters.
        effective_stiffness = self.stiffness
        effective_plasticity = self.plasticity
        effective_yield_ratio = self.yield_ratio

        if self.secondary:
            # Use the minimum age of the two particles (so the stiffer behavior of the older layer dominates).
            age = min(p1.age, p2.age)
            # Increase stiffness with age.
            effective_stiffness = self.stiffness * (1 + config.STIFFNESS_AGE_FACTOR * age)
            effective_stiffness = max(config.STIFFNESS_MIN, min(effective_stiffness, config.STIFFNESS_MAX))
            
            # Decrease plasticity with age (less plastic, i.e. more elastic).
            effective_plasticity = self.plasticity * (1 - config.PLASTICITY_AGE_FACTOR * age)
            effective_plasticity = max(config.PLASTICITY_MIN, min(effective_plasticity, config.PLASTICITY_MAX))
            
            # Increase yield_ratio with age (higher yield threshold for older layers).
            effective_yield_ratio = self.yield_ratio * (1 + config.YIELD_AGE_FACTOR * age)
            effective_yield_ratio = max(config.YIELD_RATIO_MIN, min(effective_yield_ratio, config.YIELD_RATIO_MAX))
        
        # For now, we'll use effective_stiffness in the correction.
        correction_scalar = effective_stiffness * (current_length - self.rest_length)
        correction_dir = diff / current_length
        correction = correction_scalar * correction_dir
        
        inv_mass_sum = p1.inv_mass + p2.inv_mass
        if inv_mass_sum == 0:
            return
        
        p1.pos += (p1.inv_mass / inv_mass_sum) * correction
        p2.pos -= (p2.inv_mass / inv_mass_sum) * correction

    def __repr__(self):
        return (f"Spring(particle1={self.particle1}, particle2={self.particle2}, "
                f"rest_length={self.rest_length}, stiffness={self.stiffness}, "
                f"plasticity={self.plasticity}, yield_ratio={self.yield_ratio}, "
                f"secondary={self.secondary})")
