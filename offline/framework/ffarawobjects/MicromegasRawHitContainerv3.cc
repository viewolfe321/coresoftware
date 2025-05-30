#include "MicromegasRawHitContainerv3.h"
#include "MicromegasRawHitv3.h"

#include <TClonesArray.h>

#include <iostream>

static constexpr int NHITS = 100;

MicromegasRawHitContainerv3::MicromegasRawHitContainerv3()
  : MicromegasRawHitsTCArray(new TClonesArray("MicromegasRawHitv3", NHITS))
{
}

MicromegasRawHitContainerv3::~MicromegasRawHitContainerv3()
{
  MicromegasRawHitsTCArray->Clear("C");
  delete MicromegasRawHitsTCArray;
}

void MicromegasRawHitContainerv3::Reset()
{
  MicromegasRawHitsTCArray->Clear("C");
  MicromegasRawHitsTCArray->Expand(NHITS);
}

void MicromegasRawHitContainerv3::identify(std::ostream &os) const
{
  os << "MicromegasRawHitContainerv3" << std::endl;
  os << "containing " << MicromegasRawHitsTCArray->GetEntriesFast() << " Tpc hits" << std::endl;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
  MicromegasRawHit *tpchit = static_cast<MicromegasRawHit *>(MicromegasRawHitsTCArray->At(0));
  if (tpchit)
  {
    os << "for beam clock: " << std::hex << tpchit->get_bco() << std::dec << std::endl;
  }
}

int MicromegasRawHitContainerv3::isValid() const
{
  return MicromegasRawHitsTCArray->GetSize();
}

unsigned int MicromegasRawHitContainerv3::get_nhits()
{
  return MicromegasRawHitsTCArray->GetEntriesFast();
}

MicromegasRawHit *MicromegasRawHitContainerv3::AddHit()
{
  auto *newhit = new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast() + 1]) MicromegasRawHitv3();
  return newhit;
}

MicromegasRawHit *MicromegasRawHitContainerv3::AddHit(MicromegasRawHit *rawhit)
{
  if (rawhit->IsA() == MicromegasRawHitv3::Class())
  {
    // fast add with move constructor to avoid ADC data copying
    return new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast() + 1])
        MicromegasRawHitv3(std::move(*dynamic_cast<MicromegasRawHitv3 *>(rawhit)));
  }

  // slow
  std::cout << __PRETTY_FUNCTION__ << "WARNING: input hit is not of type MicromegasRawHitv3. This is slow, please avoid." << std::endl;
  return new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast() + 1]) MicromegasRawHitv3(rawhit);
}

MicromegasRawHit *MicromegasRawHitContainerv3::get_hit(unsigned int index)
{
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
  return static_cast<MicromegasRawHit *>(MicromegasRawHitsTCArray->At(index));
}
