// Refine the grid uniformly
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::globalRefine (int refCount = 1)
{
  willCoarsen=false;

  if (maxLevel()+refCount<0)
    DUNE_THROW(GridError, "Grid has only " << maxLevel() << " levels. Cannot do "
                           << " globalRefine(" << refCount << ")");

  // The leafiterator is simply successively visiting all levels from fine to coarse
  // and just checking whether the isLeaf flag is set. Using it to identify
  // elements that need refinement would produce an endless loop, as the newly
  // added elements
  // would later be identified as elements that still need refinement.
  std::size_t oldLevels =entityImps_.size();
  //conditional typename depending on dimension of grid (1 or 2)
  typedef typename std::conditional<
                      dimgrid==2, 
                      typename std::vector<tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<1, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<2, dimgrid, dimworld> >
                                          > >::reverse_iterator,
                      typename std::vector<tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<1, dimgrid, dimworld> > 
                                          > >::reverse_iterator
                   >::type LevelIterator;

  // Allocate space for the new levels. Thus the rend iterator will not get
  // invalid due to newly added levels.
  entityImps_.reserve(oldLevels+refCount);
  levelIndexSets_.reserve(oldLevels+refCount);
  LevelIterator level = entityImps_.rbegin();

  // Add tuples for the new levels.
  for (int i=0; i < refCount; ++i)
  {
    entityImps_.push_back(EntityTuple());
    levelIndexSets_.push_back(new FoamGridLevelIndexSet<const FoamGrid >());
  }

  if (refCount < 0)
  { 
    if (globalRefined+refCount<0)
    {
      for (int i=refCount; i<0; ++i)
      {
        // Mark each leaf element for coarsening.
        typedef typename Traits::template Codim<0>::LeafIterator Iterator;
        for (Iterator elem=this->leafbegin<0>(), end = this->leafend<0>();
             elem != end; ++elem)
        {
          mark(-1,*elem);
        }

        // do the adaptation
        preAdapt();
        adapt();
        postAdapt();
      }

      globalRefined=0;
      return;
    }
    else
    {
      for (int i=refCount; i<0; ++i)
      {
        delete levelIndexSets_.back();
        levelIndexSets_.pop_back();
      }

      entityImps_.resize(oldLevels+refCount);

      // To be able to create the leaf level we need to set
      // the sons of the entities of maxlevel to null
      typename std::list<FoamGridEntityImp<0, dimgrid, dimworld> >::iterator vIt
        = Dune::get<0>(entityImps_[maxLevel()]).begin();
      typename std::list<FoamGridEntityImp<0, dimgrid, dimworld> >::iterator vEndIt
        = Dune::get<0>(entityImps_[maxLevel()]).end();
      for (; vIt!=vEndIt; ++vIt)
      {
        vIt->sons_[0]=nullptr;
        if(dimgrid == 1)
        	vIt->nSons_=0;
      }

      if(dimgrid == 2) 
      {
      	typename std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >::iterator edIt
        	= Dune::get<dimgrid-1>(entityImps_[maxLevel()]).begin();
      	typename std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >::iterator edEndIt
        	= Dune::get<dimgrid-1>(entityImps_[maxLevel()]).end();
      	for (; edIt!=edEndIt; ++edIt)
      	{
        	edIt->sons_[0]=nullptr;
        	edIt->sons_[1]=nullptr;
        	edIt->nSons_=0;
      	}
      }

      typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >::iterator elIt
        = Dune::get<dimgrid>(entityImps_[maxLevel()]).begin();
      typename std::list<FoamGridEntityImp<dimgrid, dimgrid ,dimworld> >::iterator elEndIt
        = Dune::get<dimgrid>(entityImps_[maxLevel()]).end();
      for (; elIt!=elEndIt; ++elIt)
      {
        for(int sonIdx=0; sonIdx < (1<<dimgrid); ++sonIdx) 
          elIt->sons_[sonIdx]=nullptr;
        elIt->nSons_=0;
      }
    }

  }
  else if (refCount > 0) // for refCount = 0 do nothing
  {
    // sanity check whether the above asssumption really holds.
    assert(&entityImps_[oldLevels-1]==&(*level));

    // Do the actual refinement.
    // We start with the finest level
    std::size_t levelIndex;
    for (levelIndex=oldLevels-1; level!=entityImps_.rend(); ++level, --levelIndex)
    {
      typedef typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >::iterator ElementIterator;
      bool foundLeaf=false;

      for (ElementIterator element=Dune::get<dimgrid>(*level).begin(); element != Dune::get<dimgrid>(*level).end(); ++element)
        if(element->isLeaf())
        {
          foundLeaf = true;
          dverb << "refining element " << &(*element) << std::endl;
          if (element->type().isTriangle() || element->type().isLine())
            refineSimplexElement(*element, refCount);
          else
            DUNE_THROW(NotImplemented, "Refinement only supported for simplices!");
        }

      if (!foundLeaf)
        break;
    }

    for (levelIndex+=2;levelIndex!=entityImps_.size(); ++levelIndex)
      levelIndexSets_[levelIndex]->update(*this, levelIndex);
  }

  // Update the leaf indices
  leafIndexSet_.update(*this);

  globalRefined=std::max(globalRefined+refCount,0);
  postAdapt();
}


//f Book-keeping routine to be called before adaptation
template <int dimgrid, int dimworld>
bool Dune::FoamGrid<dimgrid, dimworld>::preAdapt()
{
  // Loop over all leaf entities and check whether they might be
  // coarsened. If there is one return true.
  typedef typename Traits::template Codim<0>::LeafIterator Iterator;
  int addLevels = 0;
  willCoarsen = false;

  for (Iterator elem=this->leafbegin<0>(), end = this->leafend<0>(); elem != end; ++elem)
  {
    int mark=getMark(*elem);
    addLevels=std::max(addLevels, elem->level()+mark-maxLevel());

    if (mark<0)
    {
      // If this element is marked for coarsening, but another child
      // of this element's father is marked for refinement, then we
      // need to reset the marker to doNothing
      bool otherChildRefined=false;
      FoamGridEntityImp<dimgrid, dimgrid, dimworld>& father = *this->getRealImplementation(*elem).target_->father_;
      typedef typename array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid >::iterator ChildrenIter;
      for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
        otherChildRefined = otherChildRefined ||
                            (*child)->markState_==FoamGridEntityImp<dimgrid, dimgrid, dimworld>::REFINE;

      if (otherChildRefined)
      {
        for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
          if ((*child)->markState_==FoamGridEntityImp<dimgrid, dimgrid, dimworld>::COARSEN)
            (*child)->markState_=FoamGridEntityImp<dimgrid, dimgrid, dimworld>::DO_NOTHING;
      }
      else
        willCoarsen = willCoarsen || mark<0;
    }
  }

  if (addLevels)
  {
    entityImps_.push_back(EntityTuple());
    levelIndexSets_.push_back(new FoamGridLevelIndexSet<const FoamGrid >());
  }

  return willCoarsen;
}


// Triggers the grid refinement process
template <int dimgrid, int dimworld>
bool Dune::FoamGrid<dimgrid, dimworld>::adapt()
{
  std::set<std::size_t> levelsChanged;
  bool haveRefined=false;

  // Loop over all leaf elements and refine/coarsen those that marked for it.
  typedef typename Traits::template Codim<0>::LeafIterator Iterator;
  for (Iterator elem=this->leafbegin<0>(), end = this->leafend<0>();
       elem != end; ++elem)
  {
    int mark=getMark(*elem);
    if (mark>0)
    {
      // Refine simplices
      if (elem->type().isTriangle() || elem->type().isLine())
      {
        refineSimplexElement(*const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(*elem).target_), 1);
        levelsChanged.insert(elem->level()+1);
        haveRefined=true;
      }
      else
        DUNE_THROW(NotImplemented, "Refinement only supported for simplices!");
    }

    if (mark<0) // If simplex was already treated by coarsenSimplex mark is 0
    {
      // Coarsen simplices
      if (elem->type().isTriangle() || elem->type().isLine())
      {
        assert(elem->level());
        coarsenSimplexElement(*const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(*elem).target_));
        levelsChanged.insert(elem->level());
      }
      else
        DUNE_THROW(NotImplemented, "Coarsening only supported for simplices!");
    }
  }

  if (!willCoarsen)
    return haveRefined;

  typedef typename std::set<std::size_t>::const_reverse_iterator SIter;
  for (SIter level=levelsChanged.rbegin(); level!=levelsChanged.rend(); ++level)
  {
    // First delete the pointer
    erasePointersToEntities(Dune::get<dimgrid>(entityImps_[*level]));

    // Now delete the actual vertices.
    eraseVanishedEntities(Dune::get<0>(entityImps_[*level]));

    // erase vanished facets (only exist in 2D)
    if (dimgrid == 2)
    {
      typedef typename std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >::iterator FacetIter;
      for (FacetIter facet=Dune::get<dimgrid-1>(entityImps_[*level]).begin(),
           facetEnd=Dune::get<dimgrid-1>(entityImps_[*level]).end();
           facet != Dune::get<dimgrid-1>(entityImps_[*level]).end();)
      {
        bool hasSameLevelElements=false;
        typedef typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::const_iterator ElementIter;
        for (ElementIter elem=facet->elements_.begin();
             elem != facet->elements_.end(); ++elem)
          hasSameLevelElements = hasSameLevelElements || (*elem)->level()==*level;

        assert(facet->willVanish_!=hasSameLevelElements);

        if (!hasSameLevelElements)
        {
          assert(facet->willVanish_);
          // erase returns next position
          facet=Dune::get<dimgrid-1>(entityImps_[*level]).erase(facet);
        }
        else
          // increment
          ++facet;
      }
    }

    // And the elements
    eraseVanishedEntities(Dune::get<dimgrid>(entityImps_[*level]));

    if (Dune::get<0>(entityImps_[*level]).size())
    {
      assert(Dune::get<dimgrid-1>(entityImps_[*level]).size() &&
             Dune::get<dimgrid>(entityImps_[*level]).size());
      // Update the level indices.
      levelIndexSets_[*level]->update(*this, *level);
    }
    else
    {
      assert(!Dune::get<dimgrid-1>(entityImps_[*level]).size() &&
             !Dune::get<dimgrid>(entityImps_[*level]).size());
      if (static_cast<int>(*level)==maxLevel())
        entityImps_.pop_back();
    }
  }

  if (levelsChanged.size())
    // Update the leaf indices
    leafIndexSet_.update(*this);
  
  globalRefined=0;

  return haveRefined;
}



// Clean up refinement markers
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::postAdapt()
{
  willCoarsen=false;

  // Loop over all leaf entities and remove the isNew Marker.
  typedef typename Traits::template Codim<0>::LeafIterator Iterator;

  for (Iterator elem=this->leafbegin<0>(), end = this->leafend<0>(); elem != end; ++elem)
  {
    FoamGridEntityImp<dimgrid, dimgrid, dimworld>& element=*const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(*elem).target_);
    element.isNew_=false;
    assert(!element.willVanish_);
    if (element.father_)
      element.father_->markState_=FoamGridEntityImp<dimgrid, dimgrid, dimworld>::DO_NOTHING;
  }
}


// Erases pointers in father elements to vanished entities of the element
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::erasePointersToEntities(std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >& elements)
{
  typedef typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >::iterator EntityIterator;
  for(EntityIterator element=elements.begin();
      element != elements.end(); ++element)
  {
    if(element->willVanish_)
    {
      FoamGridEntityImp<dimgrid, dimgrid, dimworld>& father=*element->father_;

      for (unsigned int i=0; i<father.nSons_; i++)
        father.sons_[i]=nullptr;
      for (unsigned int i=0; i<father.corners(); i++)
        if (father.vertex_[i]->sons_[0]!=nullptr)
          if (father.vertex_[i]->sons_[0]->willVanish_)
            father.vertex_[i]->sons_[0]=nullptr;
      for (unsigned int i=0; i<father.corners(); i++)
        if (father.facet_[i]->sons_[0]!=nullptr
            && father.facet_[i]->sons_[0]->willVanish_)
          for (unsigned int j=0; j<dimgrid; j++)
          {
            assert(father.facet_[i]->sons_[j]!=nullptr);
            assert(father.facet_[i]->sons_[j]->willVanish_);
            father.facet_[i]->sons_[j]=nullptr;
          }
    }
  }
}

// Erase Entities from memory that vanished due to coarsening.
template <int dimgrid, int dimworld>
template <int i>
void Dune::FoamGrid<dimgrid, dimworld>::eraseVanishedEntities(std::list<FoamGridEntityImp<i, dimgrid, dimworld> >& levelEntities)
{
  typedef typename std::list<FoamGridEntityImp<i, dimgrid, dimworld> >::iterator EntityIterator;
  for (EntityIterator entity=levelEntities.begin();
      entity != levelEntities.end();)
  {
    if(entity->willVanish_)
      entity=levelEntities.erase(entity);
    else
      ++entity;
  }
}

// Coarsen a simplex element
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::coarsenSimplexElement(FoamGridEntityImp<dimgrid, dimgrid, dimworld>& element)
{
  // If we coarsen an element, this means that we erase all chidren of its father
  // to prevent inconsistencies.
  const FoamGridEntityImp<dimgrid, dimgrid, dimworld>& father = *(element.father_);

  // The facets that might be erased
  std::set<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*> childFacets;

  // The vertices that might be erased
  std::set<FoamGridEntityImp<0, dimgrid, dimworld>*> childVertices;
  typedef typename array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid >::const_iterator ChildrenIter;
  for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
  {
    // Remember element for the actual deletion taking place later
    (*child)->markState_=FoamGridEntityImp<dimgrid, dimgrid, dimworld>::IS_COARSENED;
    (*child)->willVanish_=true;

    typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::iterator FacetIter;

    for (FacetIter facet=(*child)->facet_.begin(); facet != (*child)->facet_.end(); ++facet)
    {
      // Remove references to elements that will be erased
      typedef typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::iterator
                    ElementIter;
      for (ElementIter element = (*facet)->elements_.begin();
           element != (*facet)->elements_.end();
           ++element)
      {
        for (ChildrenIter child1=father.sons_.begin(); child1 !=
             father.sons_.end(); ++child1)
          if(*element==*child1)
          {
            // To prevent removing an element from a vector
            // we just overwrite the element with its father.
            // later we will check the facets and erase all
            // of them that do not have any elements on the same
            // level.
            *element=&father;
          }
      }
    }

    typedef FoamGridEntityImp<0, dimgrid, dimworld>** VertexIter;
    for (VertexIter vertex=(*child)->vertex_; vertex != (*child)->vertex_+(*child)->corners();
        ++vertex)
      childVertices.insert(*vertex);
    typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::iterator FacetIter;
    for (FacetIter facet=(*child)->facet_.begin(); facet!= (*child)->facet_.end();
        ++facet)
      childFacets.insert(*facet);
  }

  // Check whether those guys are really erased.
  // That is, we remove all vertices and facets that are part of a child element of one
  // of the neighbours of the father element
  typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::const_iterator FacetIter;

  for(FacetIter facet=father.facet_.begin(); facet != father.facet_.end(); ++facet)
  {
    typedef typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::iterator NeighborIter;
    for (NeighborIter neighbor = (*facet)->elements_.begin();
         neighbor != (*facet)->elements_.end();
         ++neighbor)
    {
      assert((*neighbor)->level()<=father.level());
      if(*neighbor != &father && (*neighbor)->level()==father.level() &&
         !(*neighbor)->isLeaf())
      {
        // This is a real neighbor element on the same level
        // Check whether one of its children is marked for coarsening
        bool coarsened=false;
        for (ChildrenIter child=(*neighbor)->sons_.begin();
             child != (*neighbor)->sons_.end(); ++child)
          if((*child)->mightVanish())
          {
            coarsened=true;
            break;
          }

        if(!coarsened)
        {
          // Remove all entities that exist in the children of the neighbor
          // from the list
          for (ChildrenIter child=(*neighbor)->sons_.begin();
               child != (*neighbor)->sons_.end(); ++child)
          {
            typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::iterator FacetIter;
            for (FacetIter facet=(*child)->facet_.begin(); facet != (*child)->facet_.end(); ++facet)
              childFacets.erase(*facet);
            typedef FoamGridEntityImp<0, dimgrid, dimworld>** VertexIter;
            for (VertexIter vertex=(*child)->vertex_;
                 vertex != (*child)->vertex_+(*child)->corners();
                 ++vertex)
              childVertices.erase(*vertex);
          }
        }
      }
    }
  }

  typedef typename std::set<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*>::iterator SFacetIter;
  for (SFacetIter e=childFacets.begin(); e!=childFacets.end(); ++e)
    (*e)->willVanish_=true;
  typedef typename std::set<FoamGridEntityImp<0, dimgrid, dimworld>*>::iterator VertexIter;
  for (VertexIter v=childVertices.begin(); v!=childVertices.end(); ++v)
   (*v)->willVanish_=true;
}

// Refine one simplex element (2D simplex)
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::refineSimplexElement(FoamGridEntityImp<2, dimgrid, dimworld>& element,
                                                    int refCount)
{
  if(refCount<0)
  {
    // We always remove all the children from the father.
    // Removing means:
    // 1. remove pointers in father element, facets and vertices
    // 2. Mark removed entities for deletion.
    DUNE_THROW(NotImplemented, "Coarsening not implemented yet");
    return;
  }


  // TODO: Currently we are assuming that only globalRefine is available
  // and therefore the next level is not yet present.
  // For real adaptivity some of the facets and vertices might already be present.
  // Therefore we need some kind of detection for this later on.

  unsigned int nextLevel=element.level()+1;


  dvverb << "Vertices " << nextLevel << ": " << Dune::get<0>(entityImps_[nextLevel]).size()
            << std::endl;
  dvverb << "Facets " << nextLevel << ": " <<Dune::get<dimgrid-1>(entityImps_[nextLevel]).size()
            << std::endl;
  dvverb << "Elements " << nextLevel << ": " << Dune::get<dimgrid>(entityImps_[nextLevel]).size()
            << std::endl;

  array<FoamGridEntityImp<0, dimgrid, dimworld>*, 3*dimgrid> nextLevelVertices;
  std::size_t vertexIndex=0;

  // create copies of the vertices of the element
  for(unsigned int c=0; c<element.corners(); ++c)
  {
    dverb<<"Processing vertex "<<element.vertex_[c]<<std::endl;
    if(element.vertex_[c]->sons_[0]==nullptr){
      // Not refined yet
      Dune::get<0>(entityImps_[nextLevel])
      .push_back(FoamGridEntityImp<0, dimgrid, dimworld>(nextLevel,
                                               element.vertex_[c]->pos_,
                                               element.vertex_[c]->id_));
      FoamGridEntityImp<0, dimgrid, dimworld>& newVertex =
      Dune::get<0>(entityImps_[nextLevel]).back();
      element.vertex_[c]->sons_[0]=&newVertex;
    }
    check_for_duplicates(nextLevelVertices, element.vertex_[c]->sons_[0], vertexIndex);
    nextLevelVertices[vertexIndex++]=element.vertex_[c]->sons_[0];
  }

  // create new vertices from facet-midpoints together with the new facets that
  // have a father
  typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::iterator FacetIterator;

  array<FoamGridEntityImp<1, dimgrid, dimworld>*, 9> nextLevelFacets;
  std::size_t facetIndex=0;
  const Dune::ReferenceElement<double, dimgrid>& refElement
    = Dune::ReferenceElements<double, dimgrid>::general(element.type());

  // I am just to dumb for a general facet to vertice mapping.
  // Therefore we just store it here
  array<std::pair<unsigned int,unsigned int>,3 > facetVertexMapping;
  facetVertexMapping[0]=std::make_pair(0,1);
  facetVertexMapping[1]=std::make_pair(2,0);
  facetVertexMapping[2]=std::make_pair(1,2);

  for(FacetIterator facet=element.facet_.begin(); facet != element.facet_.end(); ++facet)
  {
    typedef FoamGridEntityImp<0, dimgrid, dimworld> FoamGridVertex;
    const FoamGridVertex* v0 = element.vertex_[refElement.subEntity(facetIndex/2, 1, 0, 2)];
    const FoamGridVertex* v1 = element.vertex_[refElement.subEntity(facetIndex/2, 1, 1, 2)];

    if(!(*facet)->nSons_)
    {
      // Not refined yet
      // Compute facet midpoint
      FieldVector<double, dimworld> midPoint;
      for(int dim=0; dim<dimworld;++dim)
        midPoint[dim]=((*facet)->vertex_[0]->pos_[dim]
        + (*facet)->vertex_[1]->pos_[dim]) /2.0;

      //create midpoint
      Dune::get<0>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<0, dimgrid, dimworld>(nextLevel, midPoint,
                                               freeIdCounter_[0]++));
      FoamGridEntityImp<0, dimgrid, dimworld>& midVertex =
        Dune::get<0>(entityImps_[nextLevel]).back();
      check_for_duplicates(nextLevelVertices, &midVertex, vertexIndex);
      nextLevelVertices[vertexIndex++]=&midVertex;

      // sanity check for DUNE numbering
      dvverb<<"facet "<<facetIndex/2<<": "<<"("<<(*facet)->vertex_[0]->sons_[0]<<","<<(*facet)->vertex_[1]->sons_[0]<<") with father ("<<(*facet)->vertex_[0]<<","<<(*facet)->vertex_[1]<<")"<<std::endl;
      assert(v0->sons_[0]!=nullptr);
      assert(v1->sons_[0]!=nullptr);
      assert(v0->sons_[0] == nextLevelVertices[facetVertexMapping[facetIndex/2].first] ||
        v0->sons_[0] == nextLevelVertices[facetVertexMapping[facetIndex/2].second]);

      assert(v1->sons_[0] == nextLevelVertices[facetVertexMapping[facetIndex/2].first] ||
        v1->sons_[0] == nextLevelVertices[facetVertexMapping[facetIndex/2].second]);

      // create the facets and publish them in the father
      Dune::get<1>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(nextLevelVertices[facetVertexMapping[facetIndex/2].first], &midVertex,
                                                 nextLevel, freeIdCounter_[1]++, *facet));
      (*facet)->sons_[0] = &Dune::get<1>(entityImps_[nextLevel]).back();
      ++((*facet)->nSons_);
      // Inherit the boundaryId and SegmentIndex (are treated as the same for now)
      (*facet)->sons_[0]->boundaryId_=(*facet)->boundaryId_;
      (*facet)->sons_[0]->boundarySegmentIndex_=(*facet)->boundarySegmentIndex_;
      nextLevelFacets[facetIndex++]= (*facet)->sons_[0];

      // Initialize the elements_ vector of the new facet
      // with that of the father. Later we will overwrite it
      // with the correct values.
      (*facet)->sons_[0]->elements_=(*facet)->elements_;

      assert((*facet)->vertex_[1]->sons_[0]!=nullptr);
      Dune::get<1>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(&midVertex, nextLevelVertices[facetVertexMapping[facetIndex/2].second],
                                                 nextLevel, freeIdCounter_[1]++, *facet));
      (*facet)->sons_[1] = &Dune::get<1>(entityImps_[nextLevel]).back();
      ++((*facet)->nSons_);
      // Inherit the boundaryId and SegmentIndex (are treated as the same for now)
      (*facet)->sons_[1]->boundaryId_=(*facet)->boundaryId_;
      (*facet)->sons_[1]->boundarySegmentIndex_=(*facet)->boundarySegmentIndex_;
      nextLevelFacets[facetIndex++]= (*facet)->sons_[1];

      // Initialize the elements_ vector of the new facet
      // with that of the father. Later we will overwrite it
      // with the correct values.
      (*facet)->sons_[1]->elements_=(*facet)->elements_;

      (*facet)->nSons_=2;
    } else {
      // Facets do already exist. Just add its sons to nextLevelFacets
      // but make sure that the one containing vertex facetIndex comes first
      if((*facet)->sons_[0]->vertex_[0]->id_ ==
        nextLevelVertices[facetVertexMapping[facetIndex/2].first]->id_ ||
        (*facet)->sons_[0]->vertex_[1]->id_ ==
        nextLevelVertices[facetVertexMapping[facetIndex/2].first]->id_)
      {
        nextLevelFacets[facetIndex++]=(*facet)->sons_[0];
      nextLevelFacets[facetIndex++]=(*facet)->sons_[1];
      }else{
        nextLevelFacets[facetIndex++]=(*facet)->sons_[1];
        nextLevelFacets[facetIndex++]=(*facet)->sons_[0];
      }
      if((*facet)->sons_[0]->vertex_[0]->id_!=(*facet)->vertex_[0]->id_ &&
         (*facet)->sons_[0]->vertex_[0]->id_!=(*facet)->vertex_[1]->id_)
      {
        //vertex 0 is the midpoint
        check_for_duplicates(nextLevelVertices, (*facet)->sons_[0]->vertex_[0], vertexIndex);
        nextLevelVertices[vertexIndex++]=const_cast<FoamGridEntityImp<0, dimgrid, dimworld>*>((*facet)->sons_[0]->vertex_[0]);
      } else {
        check_for_duplicates(nextLevelVertices, (*facet)->sons_[0]->vertex_[1], vertexIndex);
        nextLevelVertices[vertexIndex++]=const_cast<FoamGridEntityImp<0, dimgrid, dimworld>*>((*facet)->sons_[0]->vertex_[1]);
      }
    }
  }
  assert(facetIndex==6);
  // Create the facets that lie within the father element
  // first the one that lies opposite to the vertex 0 in the father
  Dune::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(nextLevelVertices[3],
                                             nextLevelVertices[4], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelFacets[facetIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();

  // the one opposite to father vertex 1
  Dune::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(nextLevelVertices[3],
                                             nextLevelVertices[5], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelFacets[facetIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();

  // and the one opposite to father vertex 2
  Dune::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(nextLevelVertices[4],
                                             nextLevelVertices[5], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelFacets[facetIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();

  assert(facetIndex==nextLevelFacets.size());
  assert(vertexIndex==nextLevelVertices.size());

  array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 4> nextLevelElements;
  // create the new triangles that lie in the corners
  // First the one that contains vertex 0 of the father.
  Dune::get<2>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[2]++));

  FoamGridEntityImp<dimgrid, dimgrid, dimworld>* newElement = &(Dune::get<dimgrid>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->facet_[0]=nextLevelFacets[0];
  newElement->facet_[1]=nextLevelFacets[3];
  newElement->facet_[2]=nextLevelFacets[6];
  newElement->vertex_[0]=nextLevelVertices[0];
  newElement->vertex_[1]=nextLevelVertices[3];
  newElement->vertex_[2]=nextLevelVertices[4];
  newElement->refinementIndex_=0;
  nextLevelElements[0]=newElement;
  element.sons_[0]=newElement;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;

  // Next the one that contains vertex 1 of the father.
  Dune::get<2>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(Dune::get<dimgrid>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->facet_[0]=nextLevelFacets[4];
  newElement->facet_[1]=nextLevelFacets[1];
  newElement->facet_[2]=nextLevelFacets[7];
  newElement->vertex_[0]=nextLevelVertices[1];
  newElement->vertex_[1]=nextLevelVertices[5];
  newElement->vertex_[2]=nextLevelVertices[3];
  newElement->refinementIndex_=1;
  nextLevelElements[1]=newElement;
  element.sons_[1]=newElement;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


  // Last the one that contains vertex 2 of the father.
  Dune::get<dimgrid>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid ,dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->facet_[0]=nextLevelFacets[2];
  newElement->facet_[1]=nextLevelFacets[5];
  newElement->facet_[2]=nextLevelFacets[8];
  newElement->vertex_[0]=nextLevelVertices[2];
  newElement->vertex_[1]=nextLevelVertices[4];
  newElement->vertex_[2]=nextLevelVertices[5];
  newElement->refinementIndex_=2;
  nextLevelElements[2]=newElement;
  element.sons_[2]=newElement;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


  // create the triangle in the center
  Dune::get<dimgrid>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->facet_[0]=nextLevelFacets[7];
  newElement->facet_[1]=nextLevelFacets[6];
  newElement->facet_[2]=nextLevelFacets[8];
  newElement->vertex_[0]=nextLevelVertices[3];
  newElement->vertex_[1]=nextLevelVertices[5];
  newElement->vertex_[2]=nextLevelVertices[4];
  newElement->refinementIndex_=3;
  nextLevelElements[3]=newElement;
  element.sons_[3]=newElement;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


  // Now that all the triangle are created, we can update the elements attached
  // to the facets.
  // The new (inside) neighbors of the facets lying on facets of the father element.
  std::size_t neighbors[6] = {0, 1, 2, 0, 1, 2};
  dvverb<<" element "<<&element<<std::endl;
  for(std::size_t i=0; i<6; ++i){
    // Overwrite the father element by the newly created elements.
    dvverb<<" neighbour "<<i<<": ";
    overwriteFineLevelNeighbours(*nextLevelFacets[i], nextLevelElements[neighbors[i]],
                                 &element);
  }

  // Update the neighbours of the inner facets
  nextLevelFacets[6]->elements_.push_back(nextLevelElements[0]);
  nextLevelFacets[6]->elements_.push_back(nextLevelElements[3]);
  nextLevelFacets[7]->elements_.push_back(nextLevelElements[3]);
  nextLevelFacets[7]->elements_.push_back(nextLevelElements[1]);
  nextLevelFacets[8]->elements_.push_back(nextLevelElements[3]);
  nextLevelFacets[8]->elements_.push_back(nextLevelElements[2]);

#ifndef NDEBUG
  for(std::size_t vidx=0; vidx<4; ++vidx)
  {
    dvverb<<std::endl<<" Element "<<nextLevelElements[vidx]<<": ";
    for(FacetIterator facet=nextLevelElements[vidx]->facet_.begin();
        facet != nextLevelElements[vidx]->facet_.end(); ++facet)
    {
      dvverb<<std::endl<<"   facet "<<*facet<<": ";
      typedef typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::iterator
      ElementIterator;
      bool selfFound=false;
      for(ElementIterator elem=(*facet)->elements_.begin();
          elem != (*facet)->elements_.end();
          ++elem)
      {
        dvverb << *elem<<" ";
        if(*elem == nextLevelElements[vidx])
          selfFound=true;
      }
      assert(selfFound);
    }
  }

#endif
  element.nSons_=4;

  if((refCount--)>1)
  {
    int i=0;
    typedef typename array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid>::iterator ElementIterator;
    for(ElementIterator elem=nextLevelElements.begin();
        elem != nextLevelElements.end(); ++elem)
    {
      dvverb<<std::endl<<"Refinining "<<(*elem)<<" (son of"<<&element<<") refCount="<<refCount<<" child="<<i++<<std::endl;
      refineSimplexElement(**elem, refCount);
    }
  }

  dvverb << "end refineSimplex" << std::endl;
  dvverb << "Vertices " << nextLevel << ": " << Dune::get<0>(entityImps_[nextLevel]).size()
            << std::endl;
  dvverb << "Facets " << nextLevel << ": " << Dune::get<1>(entityImps_[nextLevel]).size()
            << std::endl;
  dvverb << "Elements " << nextLevel << ": " << Dune::get<dimgrid>(entityImps_[nextLevel]).size()
            << std::endl;
}

// Refine one simplex element (1D simplex element)
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::refineSimplexElement(FoamGridEntityImp<1, dimgrid, dimworld>& element,
                                                    int refCount)
{
  if(refCount<0)
  {
    // We always remove all the children from the father.
    // Removing means:
    // 1. remove pointers in father element, facets and vertices
    // 2. Mark removed entities for deletion.
    DUNE_THROW(NotImplemented, "Coarsening not implemented yet");
    return;
  }


  // TODO: Currently we are assuming that only globalRefine is available
  // and therefore the next level is not yet present.
  // For real adaptivity some of the facets and vertices might already be present.
  // Therefore we need some kind of detection for this later on.

  unsigned int nextLevel=element.level()+1;


  dvverb << "Vertices " << nextLevel << ": " << Dune::get<0>(entityImps_[nextLevel]).size() << std::endl;
  dvverb << "Elements " << nextLevel << ": " << Dune::get<dimgrid>(entityImps_[nextLevel]).size() << std::endl;

  array<FoamGridEntityImp<0, dimgrid, dimworld>*, 3*dimgrid> nextLevelVertices;
  std::size_t vertexIndex=0;

  // create copies of the vertices of the element
  for(unsigned int c=0; c<element.corners(); ++c)
  {
    dverb<<"Processing vertex "<<element.vertex_[c]<<std::endl;
    if(element.vertex_[c]->sons_[0]==nullptr){
      // Not refined yet
      Dune::get<0>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<0, dimgrid, dimworld>(nextLevel,
                                               element.vertex_[c]->pos_,
                                               element.vertex_[c]->id_));
      FoamGridEntityImp<0, dimgrid, dimworld>& newVertex = 
        Dune::get<0>(entityImps_[nextLevel]).back();

      element.vertex_[c]->sons_[0] = &newVertex;
      // Set nSons_ in analogy to 2D grid facets
      element.vertex_[c]->nSons_++;
      assert(element.vertex_[c]->nSons_==1); // Vertex can't have more than one son
      // Inherit the boundaryId_ and the boundarySegmentIndex
      element.vertex_[c]->sons_[0]->boundaryId_= element.vertex_[c]->boundaryId_;
      element.vertex_[c]->sons_[0]->boundarySegmentIndex_= element.vertex_[c]->boundarySegmentIndex_;

      // Initialize the elements_ vector if the new facet with that of the father. Later,
      // we will overwrite it with the correct values
      element.vertex_[c]->sons_[0]->elements_= element.vertex_[c]->elements_;
    }
    //add vertex to nextLevelVertices  
    check_for_duplicates(nextLevelVertices, element.vertex_[c]->sons_[0], vertexIndex);
    nextLevelVertices[vertexIndex++]=element.vertex_[c]->sons_[0];
  }
  assert(vertexIndex==2);

  // Create the facet/vertex which lies within the father element
  // Compute element midpoint
  typedef FoamGridEntityImp<0, dimgrid, dimworld> FoamGridVertex;
  const FoamGridVertex* v0 = element.vertex_[0];
  const FoamGridVertex* v1 = element.vertex_[1];
  FieldVector<double, dimworld> midPoint;
  for(int dim=0; dim<dimworld;++dim)
         midPoint[dim]=(v0->pos_[dim] + v1->pos_[dim])*0.5;

  // Create element midpoint
  Dune::get<0>(entityImps_[nextLevel]).push_back(FoamGridEntityImp<0, dimgrid, dimworld>(nextLevel, midPoint, freeIdCounter_[0]++));
  FoamGridEntityImp<0, dimgrid, dimworld>& midVertex = Dune::get<0>(entityImps_[nextLevel]).back();
  check_for_duplicates(nextLevelVertices, &midVertex, vertexIndex);
  nextLevelVertices[vertexIndex++]=&midVertex;

  assert(vertexIndex==nextLevelVertices.size()); //==3

  // Create next level elements
  array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid > nextLevelElements;
  // create the elements
  // First the one that contains vertex 0 of the father.
  Dune::get<dimgrid>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1, 1, dimworld>(nextLevel, freeIdCounter_[dimgrid]++));

  FoamGridEntityImp<dimgrid, dimgrid, dimworld>* newElement = &(Dune::get<dimgrid>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->vertex_[0]=nextLevelVertices[0];
  newElement->vertex_[1]=nextLevelVertices[2];
  newElement->facet_[0]=nextLevelVertices[0]; // are equal to vertices but are 
  newElement->facet_[1]=nextLevelVertices[2]; // set for consistent interface
  newElement->refinementIndex_=0;
  nextLevelElements[0]=newElement;
  element.sons_[0]=newElement;
  element.nSons_++;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;

  // Next the one that contains vertex 1 of the father.
  Dune::get<dimgrid>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[dimgrid]++));
  newElement = &(Dune::get<dimgrid>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->vertex_[0]=nextLevelVertices[2];
  newElement->vertex_[1]=nextLevelVertices[1];
  newElement->facet_[0]=nextLevelVertices[2]; // are equal to vertices but are 
  newElement->facet_[1]=nextLevelVertices[1]; // set for consistent interface
  newElement->refinementIndex_=1;
  nextLevelElements[1]=newElement;
  element.sons_[1]=newElement;
  element.nSons_++;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;

  assert(element.nSons_== 1<<dimgrid); //==2

  // Now that all the elements are created, we can update the elements attached
  // to the facets.
  // The new (inside) neighbors of the facets lying on facets of the father element.
  std::size_t neighbors[2] = {0, 1};
  dvverb<<" element "<<&element<<std::endl;
  for(std::size_t i=0; i<2; ++i)
  {
    // Overwrite the father element by the newly created elements.
    dvverb<<" neighbour "<<i<<": ";
    overwriteFineLevelNeighbours(*nextLevelVertices[i], nextLevelElements[neighbors[i]],
                                    &element);
  }
  // Update the neighbours of the inner vertex
  nextLevelVertices[2]->elements_.push_back(nextLevelElements[0]);
  nextLevelVertices[2]->elements_.push_back(nextLevelElements[1]);

  if((refCount--)>1)
  {
    int i=0;
    typedef typename array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid >::iterator ElementIterator;
    for(ElementIterator elem=nextLevelElements.begin();
        elem != nextLevelElements.end(); ++elem)
    {
      dvverb<<std::endl<<"Refinining "<<(*elem)<<" (son of"<<&element<<") refCount="<<refCount<<" child="<<i++<<std::endl;
      refineSimplexElement(**elem, refCount);
    }
  }

  dvverb << "end refineSimplex" << std::endl;
  dvverb << "Vertices " << nextLevel << ": " << Dune::get<0>(entityImps_[nextLevel]).size() << std::endl;
  dvverb << "Elements " << nextLevel << ": " << Dune::get<dimgrid>(entityImps_[nextLevel]).size() << std::endl;
}

// Overwrites the neighbours of this and descendant facets
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::overwriteFineLevelNeighbours(FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>& facet,
                                                            FoamGridEntityImp<dimgrid, dimgrid, dimworld>* son,
                                                            FoamGridEntityImp<dimgrid, dimgrid, dimworld>* father)
{
  typedef typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::iterator ElementIterator;
#ifndef NDEBUG
  bool fatherFound=false;
#endif
  for(ElementIterator elem=facet.elements_.begin();
      elem != facet.elements_.end();
      ++elem)
  {
    dvverb << *elem<<" ";
    // father is replaced by the son in the elements_ vector of the vertex
    if (*elem == father)
    {
#ifndef NDEBUG
      fatherFound=true;
#endif
      *elem = son;
    }
  }

#ifndef NDEBUG
  dvverb<<std::endl;
  assert(fatherFound);
#endif

  for(std::size_t i=0; i<facet.nSons_; ++i)
    overwriteFineLevelNeighbours(*facet.sons_[i], son, father);
}


// Recompute the grid indices after the grid has changed
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::setIndices()
{
  // //////////////////////////////////////////
  //   Create the index sets
  // //////////////////////////////////////////
  for (int i=levelIndexSets_.size(); i<=maxLevel(); i++) {
    FoamGridLevelIndexSet<const FoamGrid >* p
      = new FoamGridLevelIndexSet<const FoamGrid >();
    levelIndexSets_.push_back(p);
  }

  for (int i=0; i<=maxLevel(); i++)
    if (levelIndexSets_[i])
      levelIndexSets_[i]->update(*this, i);

  // Update the leaf indices
  leafIndexSet_.update(*this);
}
