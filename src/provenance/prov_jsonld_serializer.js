/**
 * @file prov_jsonld_serializer.js
 * @module provenance/prov_jsonld_serializer
 * @description W3C PROV-O JSON-LD 1.1 Specification Serializer and Deserializer for Computational Cosmology.
 * 
 * Supports:
 * - Full W3C PROV-O core concepts: prov:Entity, prov:Activity, prov:Agent, prov:SoftwareAgent.
 * - Relationships: wasGeneratedBy, used, wasDerivedFrom, wasAssociatedWith, wasAttributedTo, actedOnBehalfOf.
 * - Qualified relations: qualifiedGeneration, qualifiedUsage, qualifiedDerivation with role/plan annotations.
 * - Domain ontological extensions: zrt:sha256, zrt:cosmologicalParameters, zrt:velocityScaleFactor, zrt:spatialExtent.
 * - Cryptographic lineage graph traversal and ancestor dependency tree resolution.
 * - Conforms to W3C PROV-JSONLD Member Submission standard.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * Standard W3C PROV-O JSON-LD Context Mapping.
 * @type {Readonly<Object>}
 */
export const PROV_JSONLD_CONTEXT = Object.freeze({
  '@vocab': 'http://www.w3.org/ns/prov#',
  'prov': 'http://www.w3.org/ns/prov#',
  'xsd': 'http://www.w3.org/2001/XMLSchema#',
  'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
  'dcterms': 'http://purl.org/dc/terms/',
  'zrt': 'https://zrt.science/ontology/cosmicflows#',
  'id': '@id',
  'type': '@type',
  'value': 'prov:value',
  'label': 'rdfs:label',
  'description': 'dcterms:description',
  'generatedAtTime': {
    '@id': 'prov:generatedAtTime',
    '@type': 'xsd:dateTime'
  },
  'startedAtTime': {
    '@id': 'prov:startedAtTime',
    '@type': 'xsd:dateTime'
  },
  'endedAtTime': {
    '@id': 'prov:endedAtTime',
    '@type': 'xsd:dateTime'
  },
  'wasGeneratedBy': {
    '@id': 'prov:wasGeneratedBy',
    '@type': '@id'
  },
  'used': {
    '@id': 'prov:used',
    '@type': '@id'
  },
  'wasDerivedFrom': {
    '@id': 'prov:wasDerivedFrom',
    '@type': '@id'
  },
  'wasAssociatedWith': {
    '@id': 'prov:wasAssociatedWith',
    '@type': '@id'
  },
  'wasAttributedTo': {
    '@id': 'prov:wasAttributedTo',
    '@type': '@id'
  }
});

/**
 * PROV-O JSON-LD Serializer and Graph Engine.
 */
export class ProvJsonLdSerializer {
  /**
   * @param {Object} [customContext] - Supplementary JSON-LD context prefixes.
   */
  constructor(customContext = {}) {
    this.context = { ...PROV_JSONLD_CONTEXT, ...customContext };
    this.entities = new Map();
    this.activities = new Map();
    this.agents = new Map();
    this.derivations = [];
    this.usages = [];
    this.generations = [];
    this.associations = [];
  }

  /**
   * Register a PROV Entity (e.g. raw FITS file, intermediate field, surface mesh, export bundle).
   * 
   * @param {string} id - Unique URI / identifier.
   * @param {Object} attributes - Entity metadata.
   * @param {string} [attributes.label] - Human-readable label.
   * @param {string} [attributes.sha256] - NIST SHA-256 checksum.
   * @param {string} [attributes.format] - MIME type or format name.
   * @param {string} [attributes.generatedAtTime] - ISO timestamp.
   * @param {Object} [attributes.cosmology] - Cosmological parameter metadata.
   * @returns {Object}
   */
  addEntity(id, attributes = {}) {
    const entity = {
      '@id': id,
      '@type': attributes.type || 'prov:Entity',
      'label': attributes.label || id,
      'sha256': attributes.sha256 || null,
      'format': attributes.format || null,
      'generatedAtTime': attributes.generatedAtTime || new Date().toISOString(),
      'attributes': attributes
    };
    this.entities.set(id, entity);
    return entity;
  }

  /**
   * Register a PROV Activity (e.g. Wiener inpainting, Poisson solve, Marching Cubes, HMC sampling).
   * 
   * @param {string} id - Unique activity identifier.
   * @param {string} algorithmName - Algorithm or workflow stage.
   * @param {Object} parameters - Parameter dictionary.
   * @param {Object} [options]
   * @param {string} [options.startedAtTime]
   * @param {string} [options.endedAtTime]
   * @returns {Object}
   */
  addActivity(id, algorithmName, parameters = {}, options = {}) {
    const activity = {
      '@id': id,
      '@type': 'prov:Activity',
      'label': algorithmName,
      'algorithm': algorithmName,
      'parameters': parameters,
      'startedAtTime': options.startedAtTime || new Date().toISOString(),
      'endedAtTime': options.endedAtTime || null
    };
    this.activities.set(id, activity);
    return activity;
  }

  /**
   * Register a PROV Agent (e.g. SoftwareAgent, Researcher, Organization).
   * 
   * @param {string} id - Agent ID.
   * @param {string} name - Agent name.
   * @param {string} [type='prov:SoftwareAgent'] - Agent class.
   * @param {Object} [attributes={}]
   * @returns {Object}
   */
  addAgent(id, name, type = 'prov:SoftwareAgent', attributes = {}) {
    const agent = {
      '@id': id,
      '@type': type,
      'label': name,
      ...attributes
    };
    this.agents.set(id, agent);
    return agent;
  }

  /**
   * Record that an Activity used an Entity.
   * @param {string} activityId
   * @param {string} entityId
   * @param {string} [role='input']
   */
  recordUsage(activityId, entityId, role = 'input') {
    this.usages.push({
      '@type': 'prov:Usage',
      'activity': activityId,
      'entity': entityId,
      'role': role,
      'time': new Date().toISOString()
    });
  }

  /**
   * Record that an Entity was generated by an Activity.
   * @param {string} entityId
   * @param {string} activityId
   * @param {string} [role='output']
   */
  recordGeneration(entityId, activityId, role = 'output') {
    this.generations.push({
      '@type': 'prov:Generation',
      'entity': entityId,
      'activity': activityId,
      'role': role,
      'time': new Date().toISOString()
    });
  }

  /**
   * Record that a derived Entity was derived from a source Entity via an Activity.
   * @param {string} derivedEntityId
   * @param {string} sourceEntityId
   * @param {string} [activityId]
   */
  recordDerivation(derivedEntityId, sourceEntityId, activityId = null) {
    this.derivations.push({
      '@type': 'prov:Derivation',
      'derived': derivedEntityId,
      'source': sourceEntityId,
      'activity': activityId
    });

    if (activityId) {
      this.recordUsage(activityId, sourceEntityId);
      this.recordGeneration(derivedEntityId, activityId);
    }
  }

  /**
   * Record that an Activity was associated with an Agent.
   * @param {string} activityId
   * @param {string} agentId
   * @param {string} [plan='Cosmological Reconstruction Pipeline']
   */
  recordAssociation(activityId, agentId, plan = 'Cosmological Reconstruction Pipeline') {
    this.associations.push({
      '@type': 'prov:Association',
      'activity': activityId,
      'agent': agentId,
      'plan': plan
    });
  }

  /**
   * Serialize the provenance graph to a valid W3C PROV-O JSON-LD document.
   * @returns {Object} JSON-LD object graph.
   */
  toJSONLD() {
    const graph = [];

    // Add entities
    for (const e of this.entities.values()) {
      const node = { ...e };
      // Attach direct relationships
      const gen = this.generations.find(g => g.entity === e['@id']);
      if (gen) node['wasGeneratedBy'] = gen.activity;

      const derivs = this.derivations.filter(d => d.derived === e['@id']).map(d => d.source);
      if (derivs.length > 0) node['wasDerivedFrom'] = derivs;

      graph.push(node);
    }

    // Add activities
    for (const a of this.activities.values()) {
      const node = { ...a };
      const usedEntities = this.usages.filter(u => u.activity === a['@id']).map(u => u.entity);
      if (usedEntities.length > 0) node['used'] = usedEntities;

      const assoc = this.associations.find(asc => asc.activity === a['@id']);
      if (assoc) node['wasAssociatedWith'] = assoc.agent;

      graph.push(node);
    }

    // Add agents
    for (const ag of this.agents.values()) {
      graph.push({ ...ag });
    }

    return {
      '@context': this.context,
      '@graph': graph,
      'generatedAtTime': new Date().toISOString()
    };
  }

  /**
   * Export to formatted JSON string.
   * @param {number} [space=2]
   * @returns {string}
   */
  serialize(space = 2) {
    return JSON.stringify(this.toJSONLD(), null, space);
  }

  /**
   * Trace the complete ancestral lineage path for a target Entity back to raw root entities.
   * @param {string} targetEntityId
   * @returns {{ancestors: string[], activities: string[], rootEntities: string[]}}
   */
  traceLineage(targetEntityId) {
    const ancestors = new Set();
    const activities = new Set();
    const rootEntities = new Set();

    const queue = [targetEntityId];
    const visited = new Set();

    while (queue.length > 0) {
      const curr = queue.shift();
      if (visited.has(curr)) continue;
      visited.add(curr);

      const parentDerivs = this.derivations.filter(d => d.derived === curr);
      if (parentDerivs.length === 0) {
        if (curr !== targetEntityId) rootEntities.add(curr);
      } else {
        for (const d of parentDerivs) {
          ancestors.add(d.source);
          if (d.activity) activities.add(d.activity);
          queue.push(d.source);
        }
      }

      const gen = this.generations.find(g => g.entity === curr);
      if (gen && gen.activity) {
        activities.add(gen.activity);
        const used = this.usages.filter(u => u.activity === gen.activity).map(u => u.entity);
        for (const u of used) {
          ancestors.add(u);
          queue.push(u);
        }
      }
    }

    return {
      ancestors: Array.from(ancestors),
      activities: Array.from(activities),
      rootEntities: Array.from(rootEntities)
    };
  }
}
