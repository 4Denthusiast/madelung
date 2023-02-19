{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Vectors (
    Lattice,
    Dimension,
    simplifyLattice,
    generateLattice,
    dual
) where

import Data.List
import GHC.TypeLits
import qualified Numeric.LinearAlgebra.Data as VD
import qualified Numeric.LinearAlgebra.Static as V
import Util (nTimes)

import Debug.Trace

-- The lattice vectors stored as columns
type Lattice (n :: Nat) = V.L n n

instance (KnownNat n, KnownNat m) => Eq (V.L n m) where
    a == b = V.extract a == V.extract b

instance KnownNat n => Eq (V.R n) where
    a == b = V.unwrap a == V.unwrap b

type Dimension n = (KnownNat n, KnownNat (n - 1), 1 <= n, ((n - 1) + 1) ~ n)

-- A reasonably small (not necessarily minimal) vector that's the sum of the given vector with a lattive vector.
simplifyByLattice :: (KnownNat m, KnownNat n) => V.L m n -> V.L m 1 -> V.L m 1
simplifyByLattice l v = v - V.mul l (V.dmmap roundDown (l V.<\> v))
    where roundDown r = fromIntegral $ round (r - 0.1 * signum r) -- Bias it towards no change in cases of numerical noise

-- A reasonably small (not necessarily minimal) equivalent representation of the same lattice.
simplifyLattice :: Dimension n => Lattice n -> Lattice n
simplifyLattice l = if l' == l then l' else simplifyLattice l'
    where l' = nTimes (fromIntegral $ natVal l) simplifyFirst l
          simplifyFirst vs = let (v,vs') = V.splitCols vs in (vs' V.||| simplifyByLattice vs' v)

-- All the lattice points within the given radius (hopefully).
generateLattice :: Dimension n => Double -> Lattice n -> [V.R n]
generateLattice r l = filter (\v -> V.dot v v <= r * r) $ map sum $ sequence $ map (\v -> let r' = fromIntegral $ ceiling (r / V.norm_2 v) in map ((v *) . V.konst) [-r'..r']) $ map (V.fromList . VD.toList) $ VD.toColumns $ V.extract $ simplifyLattice $ l

dual :: KnownNat n => Lattice n -> Lattice n
dual = V.tr . V.inv
