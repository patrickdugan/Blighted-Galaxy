# M3 checkpoint: SI combat and shearing implementation

This milestone extends [BG-PHY-02](../workpapers/physics/SI_PHYSICS.md) with a physically scaled interactive bridge.

The [GDD](../design/GDD.md) remains the gameplay contract. The source fiction retains its original text in [stories](../stories/) and [workpapers](../workpapers/).

Implemented:
- Real SI positions and delays in a playable RKM intercept scene, with 2D Lorentz velocity transformations and continuous swept collision.
- An energy-conserving *calculation* for target-frame projectile kinetic energy and a clearly fictional energy-to-shield damage coupling.
- Retarded optical position markers and proper-time telemetry.
- A 10-solar-mass black-hole and 0.65-solar-mass white-dwarf encounter with physical scale, q-based tidal diagnostic, reproducible debris and explicit mass accounting at fragmentation.
- Preservation of the previous arcade executable in a separately named binary.

Still open:
- GR near compact horizons, self-gravitating fluid breakup, causal historical ephemeris, close-contact boxing on physical units, actual collision and particle rendering benchmarking, and user playtesting.

See the detailed [bridge workpaper](../workpapers/physics/SI_GAMEPLAY_BRIDGE.md).
