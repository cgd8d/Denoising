#ifndef EventContainer_hpp
#define EventContainer_hpp

/*
Containers to hold the various information related to an event.
These objects are just passive containers.

Terminology: For scintillation denoising we can denoise a full scintillation cluster (by taking some
sort of average of the charge cluster positions) or some subset of the charge clusters which comprise
a scintillation cluster -- in the extreme, we could try to extract the scintillation energy emitted
by each charge cluster independently using the different pattern of light yields due to the different
positions where the light is emitted, though in practice we'll need to do at least some reclustering.
For charge clusters, we may be able to denoise the ionization energy of individual charge clusters or
(again) may need to do a small amount of reclustering.

As a result, we define a ScintDeposit as the collection of scintillation-emitting charge clusters which
we'll denoise as a group -- anywhere from a single charge cluster to the set of all charge clusters
associated with a single scintillation cluster.  Similarly, a ChargeDeposit is the collection of
ionization-emitting charge clusters which we'll denoise as a group -- anywhere from a single charge cluster
to some reasonable group of charge clusters which we'll treat together.  All of the charge clusters
in a ScintDeposit must be linked to the same scintillation cluster and must have a fully reconstructed
position (or else we don't have any good ideas for how to evaluate the lightmap); charge clusters in a
ChargeDeposit need only have xy position, and by default we'll treat it as coming from far away.

The Summary types contain just enough information to write the results into the output file;
these must be small because we'll be holding on to all of them until the end, and there could be
(for a long run) millions.  The other objects contain information needed by denoising itself,
such as model waveforms; these may be somewhat larger, but once we're done denoising we can throw them away.


*/

// A ScintDeposit object for us means a (set of) deposits which are treated as having
// a single scintillation energy.  This can be scintillation released at a particular point
// or the combination of scintillation released at many different points.
// The abstract contains information which lets us write a result into the output file;
// the model contains information which lets us perform denoising.
struct ScintDepositSummary
{
  double fEnergy;
  std::vector<size_t> fIncludedChargeClusters; // index of charge cluster within scint.
  size_t fScintCluster;
};

// A ChargeDeposit object for us means a (set of) deposits which are treated as having
// a single ionization energy.  Typically this will be a single charge cluster,
// though two charge clusters very close together in space and time may be merged into a single
// ChargeDeposit object.
struct ChargeDepositSummary
{
  double fMagnitude;
  std::vector<size_t> fIncludedChargeClusters;
};

// EventSummary contains the information needed to write results into the output file.
// It should be quite small, as we will keep all of them for the whole run (possibly a few million
// for long runs).
struct EventSummary
{
  int fRunNo;
  int fEventNo;
  std::vector<ScintDepositSummary> fScintDeposits;
  std::vector<ChargeDepositSummary> fChargeDeposits;
};

// EventContainer has the summary information *and* the more detailed information needed for denoising.
struct EventContainer
{
  EventSummary fEventSummary;


};
#endif
