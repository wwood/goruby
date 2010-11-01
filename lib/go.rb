require 'rsruby'
require 'bio'
require 'array_pair'

module Bio
  class Go
    def initialize
      @r = RSRuby.instance
      @r.library('GO.db')
    end
    
    # Return an array of GO identifiers that are the offspring (all the descendents)
    # of the given GO term from any ontology (cellular component, biological process
    # or molecular function)
    def go_offspring(go_id)
      o = ontology_abbreviation(go_id)
      case o
        when 'MF'
        return molecular_function_offspring(go_id)
        when 'CC'
        return cellular_component_offspring(go_id)
        when 'BP'
        return biological_process_offspring(go_id)
      else
        raise Exception, "Unknown ontology abbreviation found: #{o} for go id: #{go_id}"
      end
    end
    
    # Return an array of GO identifiers that are the offspring (all the descendents)
    # of the given GO term given that it is a cellular component
    # GO term. 
    def cellular_component_offspring(go_term)
      go_get(go_term, 'GOCCOFFSPRING')
    end
    
    # Return an array of GO identifiers that are the offspring (all the descendents)
    # of the given GO term given that it is a molecular function
    # GO term.     
    def molecular_function_offspring(go_term)
      go_get(go_term, 'GOMFOFFSPRING')
    end
    
    # Return an array of GO identifiers that are the offspring (all the descendents)
    # of the given GO term given that it is a biological process
    # GO term. 
    def biological_process_offspring(go_term)
      go_get(go_term, 'GOBPOFFSPRING')
    end
    
    # Generic method for retrieving
    # e.g offspring('GO:0042717', 'GOCCCHILDREN')
    def go_get(go_term, partition)
      answers = @r.eval_R("get('#{go_term}', #{partition})")
      return [] if answers.kind_of?(Bignum) # returns this for some reason when there's no children
      return answers
    end
    
    # Given a GO ID such as GO:0048253, return the GO term that is the 
    # primary ID (GO:0050333), so that offspring functions can be used properly.
    def primary_go_id(go_id_or_synonym_id)
      # > get('GO:0048253', GOSYNONYM)
      #GOID: GO:0050333
      #Term: thiamin-triphosphatase activity
      #Ontology: MF
      #Definition: Catalysis of the reaction: thiamin triphosphate + H2O =
      #    thiamin diphosphate + phosphate.
      #Synonym: thiamine-triphosphatase activity
      #Synonym: thiamine-triphosphate phosphohydrolase activity
      #Synonym: ThTPase activity
      #Synonym: GO:0048253
      #Secondary: GO:0048253
      
      # A performance note:
      # According to some tests that I ran, finding GOID by searching GOTERM
      # is much faster than by GOSYNONYM. A
      
      begin
        # Assume it is a primary ID, as it likely will be most of the time.
        return @r.eval_R("GOID(get('#{go_id_or_synonym_id}', GOTERM))")
      rescue RException
        # if no primary is found, try to finding it by synonym. raise RException if none is found
        begin
          return @r.eval_R("GOID(get('#{go_id_or_synonym_id}', GOSYNONYM))")
        rescue RException => e
          raise RException, "#{e.message}: GO Identifier '#{go_id_or_synonym_id}' does not appear to be a primary ID nor synonym. Is the GO.db database up to date?"
        end  
      end
    end
    
    # Retrieve the string description of the given go identifier
    def term(go_id)
      @r.eval_R("Term(get('#{go_id}', GOTERM))")
    end
    
    # Retrieve the GO annotations associated with a PDB id,
    # using Bio::Fetch PDB and UniprotKB at EBI
    def cc_pdb_to_go(pdb_id)
      # retrieve the pdb file from EBI, to extract the UniprotKB Identifiers
      pdb = Bio::Fetch.new('http://www.ebi.ac.uk/cgi-bin/dbfetch').fetch('pdb', pdb_id)
      
      # parse the PDB and return the uniprot accessions (there may be >1 because of chains)
      uniprots = Bio::PDB.new(pdb).dbref.select{|s| s.database=='UNP'}.collect{|s| s.dbAccession}
      
      gos = []
      uniprots.uniq.each do |uniprot|
        u = Bio::Fetch.new('http://www.ebi.ac.uk/cgi-bin/dbfetch').fetch('uniprot', uniprot)
        
        unp = Bio::SPTR.new(u)
        
        gos.push unp.dr('GO').select{|a|
          a['Version'].match(/^C\:/)
        }.collect{ |g|
          g['Accession']
        }
      end
      
      return gos.flatten.uniq
    end
    
    # Does the subsumer subsume the subsumee? i.e. Does it include
    # the subsumee as one of its children in the GO tree?
    # 
    # For repetitively testing one GO term subsumes others, it might
    # be faster to use subsume_tester
    def subsume?(subsumer_go_id, subsumee_go_id)
      # map the subsumee to non-synonomic id
      primaree = self.primary_go_id(subsumee_go_id)
      primarer = self.primary_go_id(subsumer_go_id)
      
      # return if they are the same - the obvious case
      return true if primaree == primarer
      
      # return if subsumee is a descendent of sumsumer
      return go_offspring(primarer).include?(primaree)
    end
    
    # Return a subsume tester for a given GO term. This method is faster
    # than repeatedly calling subsume? because the list of children is cached
    def subsume_tester(subsumer_go_id, check_for_synonym=true)
      Go::SubsumeTester.new(self, subsumer_go_id, check_for_synonym)
    end
    
    # Return 'MF', 'CC' or 'BP' corresponding to the
    def ontology_abbreviation(go_id)
      @r.eval_R("Ontology(get('#{go_id}', GOTERM))")
    end
    
    # Return an array of GO ids that correspond to the parent GO terms
    # in the ontology. This isn't the most efficient this could be, because it
    # probably gets the parents for a single id multiple times. 
    def ancestors_cc(primary_go_id)
      go_get(primary_go_id, 'GOCCANCESTOR')
    end
    
    # Return an array of ancestors of the GO term or any
    # of the GO terms' children, in no particular order. This is useful
    # when wanting to know if a term has an annotation that is 
    # non-overlapping with a particular go term. For instance, 'membrane'
    # is cordial with 'nucleus', they are boths is an ancestors of 
    # 'nuclear membrane'. However, 'mitochondrion' and 'nucleus' are 
    # not cordial, since they share no common offspring.
    def cordial_cc(primary_go_id)
      # cordial can be direct ancestors of a term - then the common term
      # is this term itself
      cordial_ids = ancestors_cc(primary_go_id)
      
      # collect all ancestors of all offspring
      offspring = cellular_component_offspring(primary_go_id)
      offspring.each do |o|
        cordial_ids.push ancestors_cc(o)
        cordial_ids.push o
      end
      
      # remove the term itself and any children - they are not
      # merely cordial
      cordial_ids = cordial_ids.flatten.uniq.reject do |i|
        offspring.include?(i) or primary_go_id==i
      end
      
      # return a uniq array of cordial terms
      cordial_ids
    end
    
    # When repeatedly testing subsumtion by a certain GO term,
    # it is faster to instantiate a SubsumeTester and use
    # Bio::GO::SubsumeTester#subsume?, rather than
    # repeatedly calling Bio::GO#subsume? because SubsumeTester
    # does caching.
    class SubsumeTester
      attr_reader :subsumer_offspring, :master_go_id
      
      def initialize(go_object, subsumer_go_id, check_for_synonym=true)
        @go = go_object
        
        if check_for_synonym
          @master_go_id = @go.primary_go_id(subsumer_go_id)
        else
          @master_go_id = subsumer_go_id
        end
        @subsumer_offspring = @go.go_offspring(@master_go_id)
        @subsumer_offspring_hash = [@subsumer_offspring].flatten.to_hash
      end
      
      def subsume?(subsumer_go_id, check_for_synonym=true)
        primaree = check_for_synonym ?
        @go.primary_go_id(subsumer_go_id) :
        subsumer_go_id
        return true if @master_go_id == primaree
        @subsumer_offspring_hash.has_key?(primaree)
      end
    end
  end
end
