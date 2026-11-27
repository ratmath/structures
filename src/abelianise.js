// G/G′ as an abelian group, for a finite Group or a presentation.
import { Presentation, presentationOf } from "./presentation.js";

/** The abelianisation G/G′ of a Group or Presentation, as an AbelianGroup (Z/d₁ × … × Z^r). */
export function abelianInvariants(G) {
  return (G instanceof Presentation ? G : presentationOf(G)).abelianInvariants();
}
